#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""tcga.py: Contains a helper class for TCGA data."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union

from numpy import iterable
from pypipegraph import Job, FileGeneratingJob
from pathlib import Path
from pandas import DataFrame
import pandas as pd
import pypipegraph as ppg
import zipfile
import requests
import tempfile
import os
import json
import subprocess
from .util import wrap

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class TCGADownloader:
    """
    take care of file downloads.
    """

    def __init__(
        self, version: str = "1.6.1", client_path: Path = Path("/project/cache/gdc_client")
    ):
        self.version = version
        self.client_file = f"gdc-client_v{version}_Ubuntu_x64.zip"
        self.client_path = client_path
        self.gdc_client_command = str(client_path / "gdc-client")

    def download_client(self) -> None:
        """
        Downloads the GDC client to the client path.
        """
        url = f"https://gdc.cancer.gov/files/public/file/{self.client_file}"
        r = requests.get(url)
        f = tempfile.TemporaryFile()
        f.write(r.content)
        with zipfile.ZipFile(f, "r") as zip_ref:
            zip_ref.extractall(self.client_path)
        os.chmod(self.gdc_client_command, 0o770)

    def ensure_client(self) -> FileGeneratingJob:
        """
        Generates a FileGeneratingJob that downloads the gdc client.

        Dependency for the actual file download.

        Returns
        -------
        FileGeneratingJob
            Job to download GDC client.
        """
        self.client_path.mkdir(parents=True, exist_ok=True)
        outfile = self.client_path / "gdc-client"
        return ppg.FileGeneratingJob(outfile, self.download_client)


class TCGAProject:
    """
    have a cases, samples, files dataframe.
    create manifest file.
    set up download.
    download missing values.
    """

    FILE_ENDPT = "https://api.gdc.cancer.gov/files"
    CASES_ENDPT = "https://api.gdc.cancer.gov/cases"
    PROJECTS_ENDPT = "https://api.gdc.cancer.gov/projects"

    def __init__(
        self, file_dir: Path = None, cache_dir: Path = None, tcgaclient: TCGADownloader = None
    ):
        self.file_dir = Path("/project/incoming/tcga")  # "/machine/ffs/datasets/tcga")
        self.cache_dir = Path("/project/cache/tcga")
        self.max_records = 100000
        self.download_client = tcgaclient if tcgaclient is not None else TCGADownloader()

    def __repr__(self):
        return f"TCGA(file_dir={self.file_dir}, cache_dir={self.cache_dir})"

    def __str__(self):
        return "TCGA"

    def fetch_files(self, df_files: DataFrame) -> List[Job]:
        self.__check_files_dataframe(df_files)
        jobs = []
        for project_id, df_files_project in df_files.groupby("project"):
            for data_type, df_files_subset in df_files_project.groupby("data_type"):
                outpath = self.file_dir / project_id / data_type
                outpath.mkdir(parents=True, exist_ok=True)
                available = self.check_available_files(outpath)
                df_files_subset_missing = df_files_subset.loc[
                    df_files_subset.index.difference(available)
                ]
                if not df_files_subset_missing.empty:
                    manifest_file = outpath / "manifest.txt"
                    manifest_job = ppg.FileGeneratingJob(
                        manifest_file,
                        wrap(self.write_manifest, df_files_subset_missing, manifest_file),
                    )
                    download_job = ppg.FileGeneratingJob(
                        manifest_file.parent / "stdout.txt",
                        wrap(self.download(manifest_file, outpath)),
                    ).dependsd_on(manifest_job)
                    jobs.extend([manifest_job, download_job])
        return jobs

    def __check_files_dataframe(self, df_files: DataFrame):
        "Checks the dataframe columns"
        missing_columns = set(["data_type", "project_id"]).difference(df_files.columns)
        if missing_columns:
            raise ValueError(
                f"The fields {missing_columns} are not present in the files DataFrame. This is needed for creating the data path."
            )
        missing_columns = set(["file_id", "file_name", "md5sum", "file_size", "state"]).difference(
            df_files.columns
        )
        if missing_columns:
            raise ValueError(
                "The following are needed for the manifest bt are missing: {missing_columns}."
            )

    def check_available_files(self, data_path: Path) -> List[str]:
        """Checks data_path for already available file ids and returns them in a List. """
        return [str(filename) for filename in data_path.iterdir() if filename.is_dir()]

    def write_manifest(self, df_files: DataFrame, outfile: Path):
        "Writes a manifest from a Dataframe."
        df = df_files[["file_id", "file_name", "md5sum", "file_size", "state"]]
        df = df.rename(columns={"file_id": "id", "file_name": "filename", "file_size": "size"})
        df.to_csv(outfile, sep="\t", index=False)

    def download(self, manifest: Path, outfolder: Path):
        """
        Starts the gdc client to download the files from a manifest.

        Parameters
        ----------
        manifest : Path
            Manifest file.
        """
        stdout = manifest.parent / "stdout.txt"
        cmd = [self.download_client.gdc_client_command] + [
            "download",
            "-m",
            str(manifest.absolute()),
            "-d",
            str(outfolder.absolute()),
        ]
        try:
            with stdout.open("w") as outp:
                subprocess.check_call(cmd, stdout=outp)
        except subprocess.CalledProcessError:
            print(" ".join(cmd))
            raise ValueError()

    def make_params(
        self, filters: str = None, fields: str = None, expand: List[str] = None, size: int = None
    ) -> Dict:
        if filters is None:
            filters = self.__get_default_filter()
        elif isinstance(filters, dict):
            filters = json.dumps(filters)
        if isinstance(fields, list):
            fields = ",".join(fields)
        if isinstance(expand, list):
            expand = ",".join(expand)
        if size is None:
            size = self.default_size
        params = {
            "filters": filters,
            "fields": fields,
            "format": "JSON",
            "size": size,
            "pretty": True,
            "expand": expand,
        }
        return params

    def __get_default_filter(self) -> str:
        filter_params = {
            "op": "=",
            "content": {"field": "cases.project.project_id", "value": [self.project_id]},
        }
        return json.dumps(filter_params)

    def __default_fields_files(self):
        fields = [
            "file_id",
            "access",
            "data_type",
            "id",
            "data_format",
            "access",
            "file_name",
            "data_category",
            "type",
            "file_size",
            "md5sum",
            "experimental_strategy",
            "cases.project.project_id",
            "cases.case_id",
            "cases.samples.sample_ids",
        ]
        return fields

    def get_files_default(self):
        fields = self.__default_fields_files()
        df = self.get_cases_default()
        return self.assert_unique_ids(df)

    def assert_unique_ids(df: DataFrame) -> DataFrame:
        assert df["id"].is_unique
        df = df.set_index("id")
        return df

    def get_cases(
        self, filters: str = None, fields: str = None, expand: List[str] = None, size: int = None
    ) -> DataFrame:
        """get a df_cases from endpoints."""
        params = self.make_params(filters, fields, expand, size)
        df = self.get_dataframe_from_endpoint(TCGAProject.CASES_ENDPT, params)
        return df

    def get_samples(self):
        """get a df_samples from endpoints."""
        pass

    def get_files(self):
        """get files df"""
        pass

    def fetch_json_from_endpoint(self, endpoint: str, params: dict):
        params["filters"] = json.dumps(params["filters"])
        response = requests.get(endpoint, params=params)
        result_json = response.json()
        return result_json

    def get_dataframe_from_endpoint(self, endpoint: str, params: dict) -> DataFrame:
        result_json = self.fetch_json_from_endpoint(endpoint, params)
        df = pd.read_json(json.dumps(result_json["data"]["hits"]))
        if df.empty:
            print("Paremeters:\n")
            print(json.dumps(params, indent=2))
            raise ValueError(
                "The parameter provided did not reeeetieve any entities from the endpoint ({endpoint})"
            )
        df = df.set_index("id")
        return df

    def check_available_files(self, data_path: Path) -> List[str]:
        """Checks data_path for already available file ids and returns them in a List. """
        with data_path.open("r") as inp:
            return [filename for filename in inp.iterdir()]

