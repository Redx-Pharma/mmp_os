#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module of unit test
"""

from MMP import utils
import pytest
import pandas as pd
import logging

log = logging.getLogger(__name__)


@pytest.fixture(scope="session")
def mmp_csv() -> str:
    df = pd.DataFrame(
        [
            ["c1ccccc1", "c1ccccc1C", "Y1", "Z1", "[c:][H] >> [c:][C]", "c1ccccc1"],
            ["c1ccccc1CF", "c1ccccc1N", "Y2", "Z2", "[c:][H] >> [c:][N]", "c1ccccc1"],
            ["c1ccccc1N", "c1ccccc1C", "Z2", "Z1", "[c:][N] >> [c:][C]", "c1ccccc1"],
            ["c1ccccc1Cl", "c1ccccc1C", "Y3", "Z1", "[c:][Cl] >> [c:][C]", "c1ccccc1"],
            ["c1ccccc1Cl", "c1ccccc1N", "Y3", "Z2", "[c:][Cl] >> [c:][N]", "c1ccccc1"],
            ["c1ccccc1C", "c1ccccc1Cl", "Y1", "Y3", "[c:][C] >> [c:][Cl]", "c1ccccc1"],
            ["CCCC", "CCCCC", "Y4", "Z3", "[C][H]>>[C][C]", "CCCC"],
            ["CCCC", "CCCCCC", "Y4", "Z4", "[C][H]>>[C][C][C]", "CCCC"],
            ["CCCCC", "CCCCC", "Y4", "Z4", "[C][H]>>[C]", "CCCCC"],
            ["C1CCCCC1", "C1CCCCC1C", "Y4", "Z4", "[C:][H]>>[C:][C]", "C1CCCCC1"],
        ],
        columns=["smiles 1", "smiles 2", "id 1", "id 2", "v1>>v2", "constant"],
    )

    df.to_csv("tmp.csv", index=False)

    return "tmp.csv"


def test_get_valid_smiles_all_valid():
    """
    Test the codes can extract the smiles and finds all valid smiles to be valid
    """
    inp = ["c1ccccc1", "CCCCC"]
    ret = utils.get_valid_smiles(inp)
    assert all(ent == ret[ith] for ith, ent in enumerate(inp))


def test_get_valid_smiles_some_invalid():
    """
    Test that invalid smiles can be remove automatically and the returned list is of the expected length
    """
    inp = ["c1ccccc1", "CCCCC", "ukn"]
    ret = utils.get_valid_smiles(inp)
    assert len(ret) == 2
    assert ret[0] == inp[0]
    assert ret[1] == inp[1]


def test_get_valid_smiles_order_some_invalid():
    """
    Test that invalid smiles can be remove automatically and the returned list is of the expected length and order
    """
    inp = ["c1ccccc1", "CCCCC", "ukn", "c1cocc1", "CCCNCC"]
    ret = utils.get_valid_smiles(inp)
    assert len(ret) == 4
    assert ret[0] == inp[0]
    assert ret[1] == inp[1]
    # Test should compare index 2 to index 3 as the invalid smiles is removed
    assert ret[2] == inp[3]
    assert ret[3] == inp[4]


def test_get_valid_smiles_mask_all_valid():
    """
    Test that we get the expected all smiles valid vector boolean mask
    """
    inp = ["c1ccccc1", "CCCCC"]
    ret = utils.get_valid_smiles_mask(inp)
    assert all(ent is True for ent in ret)


def test_get_valid_smiles_mask_some_invalid():
    """
    Test that we get the expected all but the last smiles valid vector boolean mask
    """
    inp = ["c1ccccc1", "CCCCC", "ukn"]
    ret = utils.get_valid_smiles_mask(inp)
    assert len(ret) == 3
    assert ret[0] is True
    assert ret[1] is True
    assert ret[2] is False


def test_get_valid_smiles_mask_order_some_invalid():
    """
    Test that invalid smiles are listed as False in the vector boolean mask
    """
    inp = ["c1ccccc1", "CCCCC", "ukn", "c1cocc1", "CCCNCC"]
    ret = utils.get_valid_smiles_mask(inp)
    assert len(ret) == 5
    assert ret[0] is True
    assert ret[1] is True
    assert ret[2] is False
    assert ret[3] is True
    assert ret[4] is True


def test_mmp_csv_to_mmp_dict(mmp_csv):
    """
    Function to test the building of an MMP set from a MMPDB csv
    """

    filename = mmp_csv
    mmp_data = utils.mmp_csv_to_mmp_dict(filename)

    # Matches expoected for keys : Y1 and Z2 respectively
    expected_mmp_matches = [["Z1", "Y3"], ["Y2", "Z1", "Y3"]]

    # NOTICE: The smiles is one longer as it includes the base (key) smiles as well
    expected_mmp_smiles = [
        ["c1ccccc1", "c1ccccc1C", "c1ccccc1Cl"],
        ["c1ccccc1N", "c1ccccc1CF", "c1ccccc1C", "c1ccccc1Cl"],
    ]
    log.critical(f"{mmp_data}")
    assert all(ent in mmp_data["ids"]["Y1"] for ent in expected_mmp_matches[0])
    assert all(ent in mmp_data["ids"]["Z2"] for ent in expected_mmp_matches[1])

    assert all(ent in mmp_data["smiles"]["Y1"] for ent in expected_mmp_smiles[0])
    assert all(ent in mmp_data["smiles"]["Z2"] for ent in expected_mmp_smiles[1])
