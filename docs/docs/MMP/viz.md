# Module MMP.viz

Module of methods to run MMP analysis, generation and prediction

??? example "View Source"
        #!/usr/bin.env python3

        # -*- coding: utf-8 -*-

        """

        Module of methods to run MMP analysis, generation and prediction

        """

        import logging

        from rdkit import Chem

        from rdkit.Chem import PandasTools

        from typing import List, Union, Optional

        from MMP import utils

        import numpy as np

        import pandas as pd



        log = logging.getLogger(__name__)



        def show_mmp_diffs(

            base_smiles: str,

            result: pd.DataFrame,

            smiles_column: str = "SMILES",

            structure_column: str = "structure",

            mcs_smarts: str = None,

            save_filename: Union[str, None] = None,

            highlight_overlap: bool = False,

            reduce_output: bool = True,

            column_regexes_to_reduce: Union[List[str], None] = None,

            key_property: Optional[str] = None,

            **kwargs

        ) -> pd.DataFrame:

            """

            Function to vizualize the MMP results

            Args:

                base_smiles (str): _description_

                results_df (pd.DataFrame): _description_

                smiles_column (str, optional): _description_. Defaults to "SMILES".

            Returns:

                pd.DataFrame: _description_

            """

            # work on a copy so there is no danger of overwriting

            results_df = result.copy()

            if reduce_output is True:

                if column_regexes_to_reduce is not None:

                    results_df = reduce_property_metrics(

                        results_df, search_regexs=column_regexes_to_reduce

                    )

                else:

                    results_df = reduce_property_metrics(results_df)

            # Add vizualization column

            PandasTools.AddMoleculeColumnToFrame(results_df, smiles_column, structure_column)

            base_mol = Chem.MolFromSmiles(base_smiles)

            # find mximum common sub-structure to align against

            if mcs_smarts is not None:

                mcs_mol = utils.get_maximum_common_substructure(

                    results_df[structure_column].values.tolist() + [base_mol], return_mol=True, **kwargs

                )

            else:

                mcs_mol = Chem.MolFromSmarts(mcs_smarts)

            for m in results_df[structure_column]:

                log.info(m)

                # align the molecule against the global mcs

                if m is not None:

                    PandasTools.AlignMol(m, Chem.MolToSmiles(mcs_mol))

                    # Highlight the differences compared to the base molecule using a molecule by molecule mcs

                    tmpmcs_res = utils.get_maximum_common_substructure(

                        [m, base_mol], return_mol=True

                    )

                    # copy molecule

                    # mtmp = Chem.Mol(m)

                    # highlight overlapping atoms with the base molecule

                    if highlight_overlap is True:

                        m.GetSubstructMatches(tmpmcs_res)

                    # highlight the differences with the base molecule

                    elif highlight_overlap is False:

                        matches = m.GetSubstructMatches(tmpmcs_res)

                        list_matches = [elt for tup in matches for elt in tup]

                        differences = [

                            atom.GetIdx()

                            for atom in m.GetAtoms()

                            if atom.GetIdx() not in list_matches

                        ]

                        log.debug(f"Differenes: {differences} Matches: {matches}")

                        m.__sssAtoms = differences

                else:

                    log.warning(f"Could not align molecule: {m}")

            if save_filename is not None:

                PandasTools.SaveXlsxFromFrame(

                    results_df.replace({np.nan: None}), save_filename, molCol=structure_column

                )  # , formats={"write_number": {"nan_inf_to_errors": True}})

                img = PandasTools.FrameToGridImage(

                    results_df,

                    column=structure_column,

                    molsPerRow=5,

                    subImgSize=(400, 400),

                    highlightAtomLists=[m.__sssAtoms for m in results_df["structure"]],

                    returnPNG=False,

                    useSVG=False,

                    legendsCol=key_property,

                )

                img.save(save_filename.split(".")[0] + ".png")

            return results_df



        def reduce_property_metrics(

            data_df: pd.DataFrame,

            search_regexs=[

                "_smarts",

                "_pseudosmiles",

                "_kurtosis",

                "_skewness",

                "_paired_t",

                "_p_value",

            ],

        ) -> pd.DataFrame:

            """

            Remove property columns that correspond to the regexes given as search_regexs

            Args:

                data_df (pd.DataFrame): Dataframe to remove columns from

                search_regexs (list, optional): Regexes to search for and remove columns that match them. Defaults to ["_smarts", "_pseudosmiles", "_kurtosis", "_skewness", "_q1", "_q3", "_paired_t", "_p_value"].

            Returns:

                pd.DataFrame: dataframe with column headers that match the regexs removed

            """

            df = data_df.copy()

            for regex in search_regexs:

                to_drop = [ent for ent in df.columns if ent.find(regex) != -1]

                log.debug(to_drop)

                df = df.drop(to_drop, axis=1)

            return df

## Variables

```python3
log
```

## Functions


### reduce_property_metrics

```python3
def reduce_property_metrics(
    data_df: pandas.core.frame.DataFrame,
    search_regexs=['_smarts', '_pseudosmiles', '_kurtosis', '_skewness', '_paired_t', '_p_value']
) -> pandas.core.frame.DataFrame
```

Remove property columns that correspond to the regexes given as search_regexs

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| data_df | pd.DataFrame | Dataframe to remove columns from | None |
| search_regexs | list | Regexes to search for and remove columns that match them. Defaults to ["_smarts", "_pseudosmiles", "_kurtosis", "_skewness", "_q1", "_q3", "_paired_t", "_p_value"]. | ["_smarts", "_pseudosmiles", "_kurtosis", "_skewness", "_q1", "_q3", "_paired_t", "_p_value"] |

**Returns:**

| Type | Description |
|---|---|
| pd.DataFrame | dataframe with column headers that match the regexs removed |

??? example "View Source"
        def reduce_property_metrics(

            data_df: pd.DataFrame,

            search_regexs=[

                "_smarts",

                "_pseudosmiles",

                "_kurtosis",

                "_skewness",

                "_paired_t",

                "_p_value",

            ],

        ) -> pd.DataFrame:

            """

            Remove property columns that correspond to the regexes given as search_regexs

            Args:

                data_df (pd.DataFrame): Dataframe to remove columns from

                search_regexs (list, optional): Regexes to search for and remove columns that match them. Defaults to ["_smarts", "_pseudosmiles", "_kurtosis", "_skewness", "_q1", "_q3", "_paired_t", "_p_value"].

            Returns:

                pd.DataFrame: dataframe with column headers that match the regexs removed

            """

            df = data_df.copy()

            for regex in search_regexs:

                to_drop = [ent for ent in df.columns if ent.find(regex) != -1]

                log.debug(to_drop)

                df = df.drop(to_drop, axis=1)

            return df


### show_mmp_diffs

```python3
def show_mmp_diffs(
    base_smiles: str,
    result: pandas.core.frame.DataFrame,
    smiles_column: str = 'SMILES',
    structure_column: str = 'structure',
    mcs_smarts: str = None,
    save_filename: Optional[str] = None,
    highlight_overlap: bool = False,
    reduce_output: bool = True,
    column_regexes_to_reduce: Optional[List[str]] = None,
    key_property: Optional[str] = None,
    **kwargs
) -> pandas.core.frame.DataFrame
```

Function to vizualize the MMP results

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| base_smiles | str | _description_ | None |
| results_df | pd.DataFrame | _description_ | None |
| smiles_column | str | _description_. Defaults to "SMILES". | "SMILES" |

**Returns:**

| Type | Description |
|---|---|
| pd.DataFrame | _description_ |

??? example "View Source"
        def show_mmp_diffs(

            base_smiles: str,

            result: pd.DataFrame,

            smiles_column: str = "SMILES",

            structure_column: str = "structure",

            mcs_smarts: str = None,

            save_filename: Union[str, None] = None,

            highlight_overlap: bool = False,

            reduce_output: bool = True,

            column_regexes_to_reduce: Union[List[str], None] = None,

            key_property: Optional[str] = None,

            **kwargs

        ) -> pd.DataFrame:

            """

            Function to vizualize the MMP results

            Args:

                base_smiles (str): _description_

                results_df (pd.DataFrame): _description_

                smiles_column (str, optional): _description_. Defaults to "SMILES".

            Returns:

                pd.DataFrame: _description_

            """

            # work on a copy so there is no danger of overwriting

            results_df = result.copy()

            if reduce_output is True:

                if column_regexes_to_reduce is not None:

                    results_df = reduce_property_metrics(

                        results_df, search_regexs=column_regexes_to_reduce

                    )

                else:

                    results_df = reduce_property_metrics(results_df)

            # Add vizualization column

            PandasTools.AddMoleculeColumnToFrame(results_df, smiles_column, structure_column)

            base_mol = Chem.MolFromSmiles(base_smiles)

            # find mximum common sub-structure to align against

            if mcs_smarts is not None:

                mcs_mol = utils.get_maximum_common_substructure(

                    results_df[structure_column].values.tolist() + [base_mol], return_mol=True, **kwargs

                )

            else:

                mcs_mol = Chem.MolFromSmarts(mcs_smarts)

            for m in results_df[structure_column]:

                log.info(m)

                # align the molecule against the global mcs

                if m is not None:

                    PandasTools.AlignMol(m, Chem.MolToSmiles(mcs_mol))

                    # Highlight the differences compared to the base molecule using a molecule by molecule mcs

                    tmpmcs_res = utils.get_maximum_common_substructure(

                        [m, base_mol], return_mol=True

                    )

                    # copy molecule

                    # mtmp = Chem.Mol(m)

                    # highlight overlapping atoms with the base molecule

                    if highlight_overlap is True:

                        m.GetSubstructMatches(tmpmcs_res)

                    # highlight the differences with the base molecule

                    elif highlight_overlap is False:

                        matches = m.GetSubstructMatches(tmpmcs_res)

                        list_matches = [elt for tup in matches for elt in tup]

                        differences = [

                            atom.GetIdx()

                            for atom in m.GetAtoms()

                            if atom.GetIdx() not in list_matches

                        ]

                        log.debug(f"Differenes: {differences} Matches: {matches}")

                        m.__sssAtoms = differences

                else:

                    log.warning(f"Could not align molecule: {m}")

            if save_filename is not None:

                PandasTools.SaveXlsxFromFrame(

                    results_df.replace({np.nan: None}), save_filename, molCol=structure_column

                )  # , formats={"write_number": {"nan_inf_to_errors": True}})

                img = PandasTools.FrameToGridImage(

                    results_df,

                    column=structure_column,

                    molsPerRow=5,

                    subImgSize=(400, 400),

                    highlightAtomLists=[m.__sssAtoms for m in results_df["structure"]],

                    returnPNG=False,

                    useSVG=False,

                    legendsCol=key_property,

                )

                img.save(save_filename.split(".")[0] + ".png")

            return results_df
