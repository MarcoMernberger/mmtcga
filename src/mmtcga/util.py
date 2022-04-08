#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""util.py: Contains some convenience functions needed at different locations."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union
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

    return function_wrapper()
