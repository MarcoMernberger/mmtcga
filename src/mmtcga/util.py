#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""util.py: Contains some convenience functions needed at different locations."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union
from pandas import DataFrame
import pandas as pd
import pypipegraph as ppg

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def wrap(func: Callable, *args, **kwargs) -> Callable:
    """
    Wraps a function that takes parameters and returns a parameterless
    wrapper function.

    This is needed to feed FileGeneratingJobs.

    Parameters
    ----------
    func : Callable
        Callable, taking a number of arguments/keyword arguments.

    Returns
    -------
    Callable
        A callable with bound arguments.
    """

    def function_wrapper():
        func(*args, **kwargs)

    return function_wrapper


def wrap_shared_multi(func: Callable, *args, **kwargs) -> Callable:
    """
    Wraps a function that takes parameters and returns a parameterless
    wrapper function.

    This is needed to feed FileGeneratingJobs.

    Parameters
    ----------
    func : Callable
        Callable, taking a number of arguments/keyword arguments.

    Returns
    -------
    Callable
        A callable with bound arguments.
    """

    def function_wrapper(*shared_args):
        func(*args, *shared_args)

    return function_wrapper


def explode_and_join(df: DataFrame, column_name: str) -> DataFrame:
    """explodes a column, unpacks the json dict within and joins it with the rest"""
    df_exp = df.reset_index().explode(column_name, ignore_index=True)
    df_cases = pd.json_normalize(df_exp[column_name])
    df_exp = df_exp.drop(columns=column_name, axis=1)
    assert df_cases.shape[0] == df_exp.shape[0]
    df_cases.index = df_exp.index.copy()
    doubles = df_exp.columns.intersection(df_cases.columns).values
    for c in doubles:
        if not all(df_cases[c] == df_exp[c]):
            raise ValueError(f"Columns not identical after exploding: {c}")
        df_cases = df_cases.drop(columns=c, axis=1)
    try:
        df_exp = df_exp.join(df_cases)
    except KeyError:
        print(df_exp.head(2))
        print(df_cases.head(2))
        raise
    df_exp = df_exp.set_index("id")
    return df_exp
