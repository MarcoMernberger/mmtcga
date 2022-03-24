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

    def __init__(self, project_id: str, file_dir: Path = None, cache_dir: Path = None):
        self.project_id = project_id
        self.file_dir = Path("/machine/ffs/datasets/tcga")  # "/project/incoming/tcga"
        self.cache_dir = Path("/project/cache/tcga") / self.project_id
        self.default_size = 1e7

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

