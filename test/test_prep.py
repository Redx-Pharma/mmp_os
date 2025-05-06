#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module of unit tests for the prep codes
"""

from MMP import prep
import pandas as pd
import pytest
import logging

log = logging.getLogger(__name__)


@pytest.fixture
def pandas_dataframe() -> pd.DataFrame:
    """
    pandas_dataframe fixture holding the default input dataframe

    Returns:
        pd.DataFrame - pandas data frame
    """
    test_file = pd.DataFrame(
        [
            ["c1ccccc1", "benzene", "Other"],
            ["CCCC", "butane", "Other"],
            ["NCc1ccccc1", "benzylamine", "other"],
            ["C", "methane", "something else"],
        ],
        columns=["smiles", "names", "another_column"],
    )
    return test_file


# Creates a temporty file in the global temp path
@pytest.fixture()
def to_csv_pandas_dataframe(tmp_path) -> str:
    """
    pandas_dataframe fixture holding the default input dataframe

    Returns:
        pd.DataFrame - pandas data frame
    """
    test_file = pd.DataFrame(
        [
            ["c1ccccc1", "benzene", "Other"],
            ["CCCC", "butane", "Other"],
            ["NCc1ccccc1", "benzylamine", "other"],
            ["C", "methane", "something else"],
        ],
        columns=["smiles", "names", "another_column"],
    )
    filename = tmp_path / "tmp.csv"
    test_file.to_csv(filename, index=False)

    return filename


@pytest.fixture
def pandas_dataframe_invalid_smile() -> pd.DataFrame:
    """
    pandas_dataframe_invalid_smile fixture holding the input dataframe with an invalid smiles string

    Returns:
        pd.DataFrame - pandas data frame with invalid smiles
    """
    test_file = pd.DataFrame(
        [
            ["c1ccccc1", "benzene", "Other"],
            ["CCCC", "butane", "Other"],
            ["NCc1ccccc1", "benzylamine", "other"],
            ["C", "methane", "something else"],
            [
                "NOTSMILES",
                "invalid smiles",
                "something not quite okay with these smiles",
            ],
        ],
        columns=["smiles", "names", "another_column"],
    )
    return test_file


@pytest.fixture
def expected_filtered_pandas_dataframe() -> pd.DataFrame:
    """
    expected_filtered_pandas_dataframe fixture holding the expected output dataframe from processing the pandas_dataframe for smiles and names columns only

    Returns:
        pd.DataFrame - pandas data frame with smiles and names only
    """
    test_file = pd.DataFrame(
        [
            ["c1ccccc1", "benzene"],
            ["CCCC", "butane"],
            ["NCc1ccccc1", "benzylamine"],
            ["C", "methane"],
        ],
        columns=["smiles", "names"],
    )
    return test_file


@pytest.fixture
def expected_filtered_keep_invalid_pandas_dataframe() -> pd.DataFrame:
    """
    expected_filtered_keep_invalid_pandas_dataframe fixture holding the expected output dataframe from processing the pandas_dataframe_invalid_smile for smiles and names columns only

    Returns:
        pd.DataFrame - pandas data frame with smiles and names only
    """
    test_file = pd.DataFrame(
        [
            ["c1ccccc1", "benzene"],
            ["CCCC", "butane"],
            ["NCc1ccccc1", "benzylamine"],
            ["C", "methane"],
            ["NOTSMILES", "invalid smiles"],
        ],
        columns=["smiles", "names"],
    )
    return test_file


def test_load_csv_data(pandas_dataframe, to_csv_pandas_dataframe):
    """
    test_load_csv_data test that the data can be loaded as expected from a local csv file

    Args:
        pd.DataFrame: dataframe to save to local csv
    """
    # pandas_dataframe.to_csv("tmp.csv", index=False)
    data_df = prep.load_csv_data(to_csv_pandas_dataframe)
    pd.testing.assert_frame_equal(pandas_dataframe, data_df)


def test_extract_smiles_and_labels(
    pandas_dataframe, expected_filtered_pandas_dataframe
):
    """
    test_extract_smiles_and_labels test that the data can be processed and the smiles and labels columns extracted

    Args:
        pd.DataFrame: dataframe to use as input
        pd.DataFrame: the expected pandas dataframe from the filtering
    """
    filtered_df = prep.extract_smiles_and_labels(
        pandas_dataframe, smiles_col="smiles", labels_col="names"
    )
    pd.testing.assert_frame_equal(filtered_df, expected_filtered_pandas_dataframe)


def test_extract_smiles_and_labels_with_invalid_smiles(
    pandas_dataframe_invalid_smile, expected_filtered_pandas_dataframe
):
    """
    test_extract_smiles_and_labels_with_invalid_smiles test that the data can be processed and the invalid smiles are removed
    Args:
        pd.DataFrame: dataframe to use as input
        pd.DataFrame: the expected pandas dataframe from the filtering
    """
    filtered_df = prep.extract_smiles_and_labels(
        pandas_dataframe_invalid_smile, smiles_col="smiles", labels_col="names"
    )
    pd.testing.assert_frame_equal(filtered_df, expected_filtered_pandas_dataframe)


def test_extract_smiles_and_labels_keep_invalid_smiles(
    pandas_dataframe_invalid_smile, expected_filtered_keep_invalid_pandas_dataframe
):
    """
    test_extract_smiles_and_labels_keep_invalid_smiles test that the data can be processed and the invalid smiles are kept is valid_only set to False

    Args:
        pd.DataFrame: dataframe to use as input
        pd.DataFrame: the expected pandas dataframe from the filtering
    """
    filtered_df = prep.extract_smiles_and_labels(
        pandas_dataframe_invalid_smile,
        smiles_col="smiles",
        labels_col="names",
        valid_only=False,
        rdkit_standardize_smiles=False,
    )
    pd.testing.assert_frame_equal(
        filtered_df, expected_filtered_keep_invalid_pandas_dataframe
    )
