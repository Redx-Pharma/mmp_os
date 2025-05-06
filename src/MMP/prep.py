#!/usr/bin.env python3
# -*- coding: utf-8 -*-

"""
Module of methods to prepare for MMP
"""

from typing import List
import logging
import pandas as pd
from rdkit import Chem
from MMP import utils

log = logging.getLogger(__name__)


def load_csv_data(f: str, sep: str = ",") -> pd.DataFrame:
    """
    Function to load a file of data
    Args:
        f: str file to load
        sep: str the column separator defaults to comma ","
    Returns:
        pd.DataFrame pandas dataframe of data
    """

    return pd.read_csv(f, sep=sep)


def extract_smiles_and_labels(
    data_df: pd.DataFrame,
    smiles_col: str = "smiles",
    labels_col: str = "labels",
    valid_only: bool = True,
    rdkit_standardize_smiles: bool = True,
) -> pd.DataFrame:
    """
    Function to extract the smiles and labels only for the smi file
    Args:
        data_df: dataframe of raw data including at minimum smiles_col and labels_col
        smiles_col: string smiles column name
        labels_col: string labels column name
        Valid_only: bool extract only valid smiles
    Returns:
        pandas dataframe
    """

    if valid_only is True:
        smi_valid_mask = utils.get_valid_smiles_mask(
            data_df[smiles_col].values.tolist()
        )
        df = data_df.loc[smi_valid_mask].copy()
    else:
        df = data_df.copy()

    if rdkit_standardize_smiles is True:
        df[smiles_col] = [
            Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in df[smiles_col].values
        ]

    return df[[smiles_col, labels_col]].copy()


def extract_properties(
    data_df: pd.DataFrame,
    properties_columns: List[str],
    smiles_col: str = "smiles",
    labels_col: str = "labels",
    valid_only: bool = True,
) -> pd.DataFrame:
    """
    Function to extract the smiles and labels only for the smi file
    Args:
        data_df: dataframe of raw data including at minimum smiles_col and labels_col
        smiles_col: string smiles column name
        labels_col: string labels column name
        Valid_only: bool extract only valid smiles
    Returns:
        pandas dataframe
    """

    if valid_only is True:
        smi_valid_mask = utils.get_valid_smiles_mask(
            data_df[smiles_col].values.tolist()
        )
        df = data_df.loc[smi_valid_mask].copy()
    else:
        df = data_df.copy()

    df = df[[labels_col] + list(properties_columns)].copy()
    df = df.rename({labels_col: "ID"}, axis=1)
    df = df.fillna("*")

    return df
