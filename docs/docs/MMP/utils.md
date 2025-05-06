# Module MMP.utils

Utilities module

??? example "View Source"
        #!/usr/bin.env python3

        # -*- coding: utf-8 -*-

        """

        Utilities module

        """

        from typing import List, Union, Optional

        import rdkit

        import pandas as pd

        from rdkit import Chem

        from rdkit.Chem import MCS

        from rdkit.Chem import rdFMCS

        import logging

        log = logging.getLogger(__name__)



        def get_valid_smiles(smis: List[str]) -> List[str]:

            """

            Function to extract only valid (RDKit parsable) smiles strings from a list of smiles strings

            Args:

                smis (List[str]): list of smiles strings to determine which are valid

            Returns:

                List[str]: valid (RDKit parsable) smiles strings

            """

            log.info(

                f"Determining smiles validity using RDKit parsing for RDKit version {rdkit.__version__}"

            )

            return [s for s in smis if Chem.MolFromSmiles(s) is not None]



        def get_valid_smiles_mask(smis: List[str]) -> List[str]:

            """

            Function to extract a boolean vector defining whether the smiles string in that row is valid (RDKit parsable) or not

            Args:

                smis (List[str]): list of smiles strings to determine which are valid

            Returns:

                List[str]: boolean vector defining if the smiles string of that row index is valid (RDKit parsable)

            """

            log.info(

                f"Determining smiles validity using RDKit parsing for RDKit version {rdkit.__version__}"

            )

            return [True if Chem.MolFromSmiles(s) is not None else False for s in smis]



        def get_maximum_common_substructure(

            mols: List[Chem.rdchem.Mol],

            return_smarts: bool = False,

            return_smiles: bool = False,

            return_mol: bool = False,

            print_results: bool = False,

            **kwargs,

        ) -> Union[str, Chem.rdchem.Mol, rdkit.Chem.rdFMCS.MCSResult]:

            """

            Function to find the maximum common substructure across a list of molecules

            Args:

                mols (List[rdkit.Chem.rdchem.Mol]): List of RDKit molecule objects

            Returns:

                Based on user request:

                    * (Default): rdkit.Chem.MCS.MCSResult The raw MCS result

                    * return_smart is True: str Smarts string

                    * return_smiles is True: str Smiles string

                    * return_mol is True: Chem.rdchem.Mol Molecule object

            Doctest:

            >>> get_maximum_common_substructure([Chem.MolFromSmiles("Cc1ccccc1"), Chem.MolFromSmiles("Nc1ccccc1"), Chem.MolFromSmiles("Oc1ccccc1"), Chem.MolFromSmiles("CCc1ccccc1")], return_smiles=True)

            'C1:C:C:C:C:C:1'

            >>> get_maximum_common_substructure([Chem.MolFromSmiles("Cc1ccccc1"), Chem.MolFromSmiles("Nc1ccccc1"), Chem.MolFromSmiles("Oc1ccccc1"), Chem.MolFromSmiles("CCc1ccccc1")], return_smarts=True)

            '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'

            """

            mcs_res = rdFMCS.FindMCS(mols, **kwargs)

            if print_results is True:

                log.info(f"Found MCS: {mcs_res.smartsString}")

            if return_smarts:

                return mcs_res.smartsString

            elif return_smiles:

                return Chem.MolToSmiles(Chem.MolFromSmarts(mcs_res.smartsString))

            elif return_mol:

                return Chem.MolFromSmarts(mcs_res.smartsString)

            else:

                return mcs_res

        # def get_maximum_common_substructure(

        #     mols: List[Chem.rdchem.Mol],

        #     return_smarts: bool = False,

        #     return_smiles: bool = False,

        #     return_mol: bool = False,

        # ) -> Union[str, Chem.rdchem.Mol, Chem.MCS.MCSResult]:

        #     """

        #     Function to find the maximum common substructure across a list of molecules

        #     Args:

        #         mols (List[rdkit.Chem.rdchem.Mol]): List of RDKit molecule objects

        #     Returns:

        #         Based on user request:

        #             * (Default): rdkit.Chem.MCS.MCSResult The raw MCS result

        #             * return_smart is True: str Smarts string

        #             * return_smiles is True: str Smiles string

        #             * return_mol is True: Chem.rdchem.Mol Molecule object

        #     """

        #     mcs_res = MCS.FindMCS(mols)

        #     if return_smarts:

        #         return mcs_res.smarts

        #     elif return_smiles:

        #         return Chem.MolToSmiles(Chem.MolFromSmarts(mcs_res.smarts))

        #     elif return_mol:

        #         return Chem.MolFromSmarts(mcs_res.smarts)

        #     else:

        #         return mcs_res



        get_mcs = get_maximum_common_substructure



        def mmp_csv_to_mmp_dict_with_duplicates(

            mmp_csv: Union[str, pd.DataFrame],

            id_column_one: Optional[str] = "id 1",

            id_column_two: Optional[str] = "id 2",

            smiles_column_one: Optional[str] = "smiles 1",

            smiles_column_two: Optional[str] = "smiles 2",

        ) -> dict:

            """

            Function to map a mmpdb index run output in csv format to a python dictionary of matched molecular pairs. Duplicates can arise from multiple matched changes between two molecules.

            Args:

                mmp_csv (Union[str, pd.DataFrame]): dataset from MMPDB in csv format

                id_column_one (Optional[str], optional): Matched pair id column 1. Defaults to "id 1".

                id_column_two (Optional[str], optional): Matched pair id column 2. Defaults to "id 2".

                smiles_column_one (Optional[str], optional): Matched pair smiles column 1. Defaults to "smiles 1".

                smiles_column_two (Optional[str], optional): Matched pair smiles column 1. Defaults to "smiles 2".

            Raises:

                TypeError: if incorrect input type is given

            Returns:

                dict: matched pairs ids and smiles. The ids are mapped key id : [list of all matching ids]. The smiles are mapped key id: [key id smiles] + [list of all matching ids smiles]

            """

            if isinstance(mmp_csv, str):

                tmp_df = pd.read_csv(mmp_csv)

            elif isinstance(mmp_csv, pd.DataFrame):

                tmp_df = mmp_csv.copy()

            else:

                raise TypeError(

                    f"mmp_csv must be a string to a csv file to load or a pandas dataframe not a {type(mmp_csv)}"

                )

            mmp_dict = {"ids": {}, "smiles": {}}

            for elt in set(tmp_df[id_column_one].to_list()):

                log.info(

                    f"{elt} occrus in the id 1 list {len([ent for ent in tmp_df['id 1'].to_list() if ent  == elt])} times"

                )

                base_smiles = tmp_df[tmp_df[id_column_one] == elt][smiles_column_one].tolist()[

                    0

                ]

                mmp_dict["ids"][elt] = [

                    tmp_df.loc[ith, id_column_two]

                    for ith, ent in enumerate(tmp_df[id_column_one].to_list())

                    if ent == elt

                ]

                mmp_dict["smiles"][elt] = [base_smiles] + [

                    tmp_df.loc[ith, smiles_column_two]

                    for ith, ent in enumerate(tmp_df[id_column_one].to_list())

                    if ent == elt

                ]

            return mmp_dict



        def mmp_csv_to_mmp_dict(

            mmp_csv: Union[str, pd.DataFrame],

            id_column_one: Optional[str] = "id 1",

            id_column_two: Optional[str] = "id 2",

            smiles_column_one: Optional[str] = "smiles 1",

            smiles_column_two: Optional[str] = "smiles 2",

        ) -> dict:

            """

            Function to map a mmpdb index run output in csv format to a python dictionary of matched molecular pairs

            Args:

                mmp_csv (Union[str, pd.DataFrame]): dataset from MMPDB in csv format

                id_column_one (Optional[str], optional): Matched pair id column 1. Defaults to "id 1".

                id_column_two (Optional[str], optional): Matched pair id column 2. Defaults to "id 2".

                smiles_column_one (Optional[str], optional): Matched pair smiles column 1. Defaults to "smiles 1".

                smiles_column_two (Optional[str], optional): Matched pair smiles column 1. Defaults to "smiles 2".

            Raises:

                TypeError: if incorrect input type is given

            Returns:

                dict: matched pairs ids and smiles. The ids are mapped key id : [list of all matching ids de-duplicated]. The smiles are mapped key id: [key id smiles] + [list of all matching ids smiles de-duplicated]

            """

            if isinstance(mmp_csv, str):

                tmp_df = pd.read_csv(mmp_csv)

            elif isinstance(mmp_csv, pd.DataFrame):

                tmp_df = mmp_csv.copy()

            else:

                raise TypeError(

                    f"mmp_csv must be a string to a csv file to load or a pandas dataframe not a {type(mmp_csv)}"

                )

            mmp_dict = {"ids": {}, "smiles": {}}

            for elt in set(tmp_df[id_column_one].to_list() + tmp_df[id_column_two].to_list()):

                log.info(

                    f"{elt} occrus in the id 1 list {len([ent for ent in tmp_df[id_column_one].to_list() if ent  == elt])} times and list 2 {len([ent for ent in tmp_df[id_column_two].to_list() if ent  == elt])} times"

                )

                try:

                    base_smiles = tmp_df[tmp_df[id_column_one] == elt][

                        smiles_column_one

                    ].tolist()[0]

                except IndexError:

                    base_smiles = tmp_df[tmp_df[id_column_two] == elt][

                        smiles_column_two

                    ].tolist()[0]

                ids = []

                smiles = []

                for ith, ent in enumerate(tmp_df[id_column_one].to_list()):

                    if (

                        ent.strip().lower() == elt.strip().lower()

                        and tmp_df.loc[ith, id_column_two].strip() not in ids

                    ):

                        ids.append(tmp_df.loc[ith, id_column_two].strip())

                        smiles.append(tmp_df.loc[ith, smiles_column_two].strip())

                for ith, ent in enumerate(tmp_df[id_column_two].to_list()):

                    if (

                        ent.strip().lower() == elt.strip().lower()

                        and tmp_df.loc[ith, id_column_one].strip() not in ids

                    ):

                        ids.append(tmp_df.loc[ith, id_column_one].strip())

                        smiles.append(tmp_df.loc[ith, smiles_column_one].strip())

                mmp_dict["ids"][elt] = ids

                mmp_dict["smiles"][elt] = [base_smiles] + smiles

            return mmp_dict

## Variables

```python3
log
```

## Functions


### get_maximum_common_substructure

```python3
def get_maximum_common_substructure(
    mols: List[rdkit.Chem.rdchem.Mol],
    return_smarts: bool = False,
    return_smiles: bool = False,
    return_mol: bool = False,
    print_results: bool = False,
    **kwargs
) -> Union[str, rdkit.Chem.rdchem.Mol, rdkit.Chem.rdFMCS.MCSResult]
```

Function to find the maximum common substructure across a list of molecules

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| mols | List[rdkit.Chem.rdchem.Mol] | List of RDKit molecule objects | None |

**Returns:**

| Type | Description |
|---|---|
| None | Based on user request:<br>* (Default): rdkit.Chem.MCS.MCSResult The raw MCS result<br>* return_smart is True: str Smarts string<br>* return_smiles is True: str Smiles string<br>* return_mol is True: Chem.rdchem.Mol Molecule object |

??? example "View Source"
        def get_maximum_common_substructure(

            mols: List[Chem.rdchem.Mol],

            return_smarts: bool = False,

            return_smiles: bool = False,

            return_mol: bool = False,

            print_results: bool = False,

            **kwargs,

        ) -> Union[str, Chem.rdchem.Mol, rdkit.Chem.rdFMCS.MCSResult]:

            """

            Function to find the maximum common substructure across a list of molecules

            Args:

                mols (List[rdkit.Chem.rdchem.Mol]): List of RDKit molecule objects

            Returns:

                Based on user request:

                    * (Default): rdkit.Chem.MCS.MCSResult The raw MCS result

                    * return_smart is True: str Smarts string

                    * return_smiles is True: str Smiles string

                    * return_mol is True: Chem.rdchem.Mol Molecule object

            Doctest:

            >>> get_maximum_common_substructure([Chem.MolFromSmiles("Cc1ccccc1"), Chem.MolFromSmiles("Nc1ccccc1"), Chem.MolFromSmiles("Oc1ccccc1"), Chem.MolFromSmiles("CCc1ccccc1")], return_smiles=True)

            'C1:C:C:C:C:C:1'

            >>> get_maximum_common_substructure([Chem.MolFromSmiles("Cc1ccccc1"), Chem.MolFromSmiles("Nc1ccccc1"), Chem.MolFromSmiles("Oc1ccccc1"), Chem.MolFromSmiles("CCc1ccccc1")], return_smarts=True)

            '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'

            """

            mcs_res = rdFMCS.FindMCS(mols, **kwargs)

            if print_results is True:

                log.info(f"Found MCS: {mcs_res.smartsString}")

            if return_smarts:

                return mcs_res.smartsString

            elif return_smiles:

                return Chem.MolToSmiles(Chem.MolFromSmarts(mcs_res.smartsString))

            elif return_mol:

                return Chem.MolFromSmarts(mcs_res.smartsString)

            else:

                return mcs_res


### get_mcs

```python3
def get_mcs(
    mols: List[rdkit.Chem.rdchem.Mol],
    return_smarts: bool = False,
    return_smiles: bool = False,
    return_mol: bool = False,
    print_results: bool = False,
    **kwargs
) -> Union[str, rdkit.Chem.rdchem.Mol, rdkit.Chem.rdFMCS.MCSResult]
```

Function to find the maximum common substructure across a list of molecules

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| mols | List[rdkit.Chem.rdchem.Mol] | List of RDKit molecule objects | None |

**Returns:**

| Type | Description |
|---|---|
| None | Based on user request:<br>* (Default): rdkit.Chem.MCS.MCSResult The raw MCS result<br>* return_smart is True: str Smarts string<br>* return_smiles is True: str Smiles string<br>* return_mol is True: Chem.rdchem.Mol Molecule object |

??? example "View Source"
        def get_maximum_common_substructure(

            mols: List[Chem.rdchem.Mol],

            return_smarts: bool = False,

            return_smiles: bool = False,

            return_mol: bool = False,

            print_results: bool = False,

            **kwargs,

        ) -> Union[str, Chem.rdchem.Mol, rdkit.Chem.rdFMCS.MCSResult]:

            """

            Function to find the maximum common substructure across a list of molecules

            Args:

                mols (List[rdkit.Chem.rdchem.Mol]): List of RDKit molecule objects

            Returns:

                Based on user request:

                    * (Default): rdkit.Chem.MCS.MCSResult The raw MCS result

                    * return_smart is True: str Smarts string

                    * return_smiles is True: str Smiles string

                    * return_mol is True: Chem.rdchem.Mol Molecule object

            Doctest:

            >>> get_maximum_common_substructure([Chem.MolFromSmiles("Cc1ccccc1"), Chem.MolFromSmiles("Nc1ccccc1"), Chem.MolFromSmiles("Oc1ccccc1"), Chem.MolFromSmiles("CCc1ccccc1")], return_smiles=True)

            'C1:C:C:C:C:C:1'

            >>> get_maximum_common_substructure([Chem.MolFromSmiles("Cc1ccccc1"), Chem.MolFromSmiles("Nc1ccccc1"), Chem.MolFromSmiles("Oc1ccccc1"), Chem.MolFromSmiles("CCc1ccccc1")], return_smarts=True)

            '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'

            """

            mcs_res = rdFMCS.FindMCS(mols, **kwargs)

            if print_results is True:

                log.info(f"Found MCS: {mcs_res.smartsString}")

            if return_smarts:

                return mcs_res.smartsString

            elif return_smiles:

                return Chem.MolToSmiles(Chem.MolFromSmarts(mcs_res.smartsString))

            elif return_mol:

                return Chem.MolFromSmarts(mcs_res.smartsString)

            else:

                return mcs_res


### get_valid_smiles

```python3
def get_valid_smiles(
    smis: List[str]
) -> List[str]
```

Function to extract only valid (RDKit parsable) smiles strings from a list of smiles strings

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| smis | List[str] | list of smiles strings to determine which are valid | None |

**Returns:**

| Type | Description |
|---|---|
| List[str] | valid (RDKit parsable) smiles strings |

??? example "View Source"
        def get_valid_smiles(smis: List[str]) -> List[str]:

            """

            Function to extract only valid (RDKit parsable) smiles strings from a list of smiles strings

            Args:

                smis (List[str]): list of smiles strings to determine which are valid

            Returns:

                List[str]: valid (RDKit parsable) smiles strings

            """

            log.info(

                f"Determining smiles validity using RDKit parsing for RDKit version {rdkit.__version__}"

            )

            return [s for s in smis if Chem.MolFromSmiles(s) is not None]


### get_valid_smiles_mask

```python3
def get_valid_smiles_mask(
    smis: List[str]
) -> List[str]
```

Function to extract a boolean vector defining whether the smiles string in that row is valid (RDKit parsable) or not

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| smis | List[str] | list of smiles strings to determine which are valid | None |

**Returns:**

| Type | Description |
|---|---|
| List[str] | boolean vector defining if the smiles string of that row index is valid (RDKit parsable) |

??? example "View Source"
        def get_valid_smiles_mask(smis: List[str]) -> List[str]:

            """

            Function to extract a boolean vector defining whether the smiles string in that row is valid (RDKit parsable) or not

            Args:

                smis (List[str]): list of smiles strings to determine which are valid

            Returns:

                List[str]: boolean vector defining if the smiles string of that row index is valid (RDKit parsable)

            """

            log.info(

                f"Determining smiles validity using RDKit parsing for RDKit version {rdkit.__version__}"

            )

            return [True if Chem.MolFromSmiles(s) is not None else False for s in smis]


### mmp_csv_to_mmp_dict

```python3
def mmp_csv_to_mmp_dict(
    mmp_csv: Union[str, pandas.core.frame.DataFrame],
    id_column_one: Optional[str] = 'id 1',
    id_column_two: Optional[str] = 'id 2',
    smiles_column_one: Optional[str] = 'smiles 1',
    smiles_column_two: Optional[str] = 'smiles 2'
) -> dict
```

Function to map a mmpdb index run output in csv format to a python dictionary of matched molecular pairs

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| mmp_csv | Union[str, pd.DataFrame] | dataset from MMPDB in csv format | None |
| id_column_one | Optional[str] | Matched pair id column 1. Defaults to "id 1". | "id 1" |
| id_column_two | Optional[str] | Matched pair id column 2. Defaults to "id 2". | "id 2" |
| smiles_column_one | Optional[str] | Matched pair smiles column 1. Defaults to "smiles 1". | "smiles 1" |
| smiles_column_two | Optional[str] | Matched pair smiles column 1. Defaults to "smiles 2". | "smiles 2" |

**Returns:**

| Type | Description |
|---|---|
| dict | matched pairs ids and smiles. The ids are mapped key id : [list of all matching ids de-duplicated]. The smiles are mapped key id: [key id smiles] + [list of all matching ids smiles de-duplicated] |

**Raises:**

| Type | Description |
|---|---|
| TypeError | if incorrect input type is given |

??? example "View Source"
        def mmp_csv_to_mmp_dict(

            mmp_csv: Union[str, pd.DataFrame],

            id_column_one: Optional[str] = "id 1",

            id_column_two: Optional[str] = "id 2",

            smiles_column_one: Optional[str] = "smiles 1",

            smiles_column_two: Optional[str] = "smiles 2",

        ) -> dict:

            """

            Function to map a mmpdb index run output in csv format to a python dictionary of matched molecular pairs

            Args:

                mmp_csv (Union[str, pd.DataFrame]): dataset from MMPDB in csv format

                id_column_one (Optional[str], optional): Matched pair id column 1. Defaults to "id 1".

                id_column_two (Optional[str], optional): Matched pair id column 2. Defaults to "id 2".

                smiles_column_one (Optional[str], optional): Matched pair smiles column 1. Defaults to "smiles 1".

                smiles_column_two (Optional[str], optional): Matched pair smiles column 1. Defaults to "smiles 2".

            Raises:

                TypeError: if incorrect input type is given

            Returns:

                dict: matched pairs ids and smiles. The ids are mapped key id : [list of all matching ids de-duplicated]. The smiles are mapped key id: [key id smiles] + [list of all matching ids smiles de-duplicated]

            """

            if isinstance(mmp_csv, str):

                tmp_df = pd.read_csv(mmp_csv)

            elif isinstance(mmp_csv, pd.DataFrame):

                tmp_df = mmp_csv.copy()

            else:

                raise TypeError(

                    f"mmp_csv must be a string to a csv file to load or a pandas dataframe not a {type(mmp_csv)}"

                )

            mmp_dict = {"ids": {}, "smiles": {}}

            for elt in set(tmp_df[id_column_one].to_list() + tmp_df[id_column_two].to_list()):

                log.info(

                    f"{elt} occrus in the id 1 list {len([ent for ent in tmp_df[id_column_one].to_list() if ent  == elt])} times and list 2 {len([ent for ent in tmp_df[id_column_two].to_list() if ent  == elt])} times"

                )

                try:

                    base_smiles = tmp_df[tmp_df[id_column_one] == elt][

                        smiles_column_one

                    ].tolist()[0]

                except IndexError:

                    base_smiles = tmp_df[tmp_df[id_column_two] == elt][

                        smiles_column_two

                    ].tolist()[0]

                ids = []

                smiles = []

                for ith, ent in enumerate(tmp_df[id_column_one].to_list()):

                    if (

                        ent.strip().lower() == elt.strip().lower()

                        and tmp_df.loc[ith, id_column_two].strip() not in ids

                    ):

                        ids.append(tmp_df.loc[ith, id_column_two].strip())

                        smiles.append(tmp_df.loc[ith, smiles_column_two].strip())

                for ith, ent in enumerate(tmp_df[id_column_two].to_list()):

                    if (

                        ent.strip().lower() == elt.strip().lower()

                        and tmp_df.loc[ith, id_column_one].strip() not in ids

                    ):

                        ids.append(tmp_df.loc[ith, id_column_one].strip())

                        smiles.append(tmp_df.loc[ith, smiles_column_one].strip())

                mmp_dict["ids"][elt] = ids

                mmp_dict["smiles"][elt] = [base_smiles] + smiles

            return mmp_dict


### mmp_csv_to_mmp_dict_with_duplicates

```python3
def mmp_csv_to_mmp_dict_with_duplicates(
    mmp_csv: Union[str, pandas.core.frame.DataFrame],
    id_column_one: Optional[str] = 'id 1',
    id_column_two: Optional[str] = 'id 2',
    smiles_column_one: Optional[str] = 'smiles 1',
    smiles_column_two: Optional[str] = 'smiles 2'
) -> dict
```

Function to map a mmpdb index run output in csv format to a python dictionary of matched molecular pairs. Duplicates can arise from multiple matched changes between two molecules.

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| mmp_csv | Union[str, pd.DataFrame] | dataset from MMPDB in csv format | None |
| id_column_one | Optional[str] | Matched pair id column 1. Defaults to "id 1". | "id 1" |
| id_column_two | Optional[str] | Matched pair id column 2. Defaults to "id 2". | "id 2" |
| smiles_column_one | Optional[str] | Matched pair smiles column 1. Defaults to "smiles 1". | "smiles 1" |
| smiles_column_two | Optional[str] | Matched pair smiles column 1. Defaults to "smiles 2". | "smiles 2" |

**Returns:**

| Type | Description |
|---|---|
| dict | matched pairs ids and smiles. The ids are mapped key id : [list of all matching ids]. The smiles are mapped key id: [key id smiles] + [list of all matching ids smiles] |

**Raises:**

| Type | Description |
|---|---|
| TypeError | if incorrect input type is given |

??? example "View Source"
        def mmp_csv_to_mmp_dict_with_duplicates(

            mmp_csv: Union[str, pd.DataFrame],

            id_column_one: Optional[str] = "id 1",

            id_column_two: Optional[str] = "id 2",

            smiles_column_one: Optional[str] = "smiles 1",

            smiles_column_two: Optional[str] = "smiles 2",

        ) -> dict:

            """

            Function to map a mmpdb index run output in csv format to a python dictionary of matched molecular pairs. Duplicates can arise from multiple matched changes between two molecules.

            Args:

                mmp_csv (Union[str, pd.DataFrame]): dataset from MMPDB in csv format

                id_column_one (Optional[str], optional): Matched pair id column 1. Defaults to "id 1".

                id_column_two (Optional[str], optional): Matched pair id column 2. Defaults to "id 2".

                smiles_column_one (Optional[str], optional): Matched pair smiles column 1. Defaults to "smiles 1".

                smiles_column_two (Optional[str], optional): Matched pair smiles column 1. Defaults to "smiles 2".

            Raises:

                TypeError: if incorrect input type is given

            Returns:

                dict: matched pairs ids and smiles. The ids are mapped key id : [list of all matching ids]. The smiles are mapped key id: [key id smiles] + [list of all matching ids smiles]

            """

            if isinstance(mmp_csv, str):

                tmp_df = pd.read_csv(mmp_csv)

            elif isinstance(mmp_csv, pd.DataFrame):

                tmp_df = mmp_csv.copy()

            else:

                raise TypeError(

                    f"mmp_csv must be a string to a csv file to load or a pandas dataframe not a {type(mmp_csv)}"

                )

            mmp_dict = {"ids": {}, "smiles": {}}

            for elt in set(tmp_df[id_column_one].to_list()):

                log.info(

                    f"{elt} occrus in the id 1 list {len([ent for ent in tmp_df['id 1'].to_list() if ent  == elt])} times"

                )

                base_smiles = tmp_df[tmp_df[id_column_one] == elt][smiles_column_one].tolist()[

                    0

                ]

                mmp_dict["ids"][elt] = [

                    tmp_df.loc[ith, id_column_two]

                    for ith, ent in enumerate(tmp_df[id_column_one].to_list())

                    if ent == elt

                ]

                mmp_dict["smiles"][elt] = [base_smiles] + [

                    tmp_df.loc[ith, smiles_column_two]

                    for ith, ent in enumerate(tmp_df[id_column_one].to_list())

                    if ent == elt

                ]

            return mmp_dict
