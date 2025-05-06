# Module MMP.runmmp

Module of methods to run MMP analysis, generation and prediction

??? example "View Source"
        #!/usr/bin.env python3

        # -*- coding: utf-8 -*-

        """

        Module of methods to run MMP analysis, generation and prediction

        """

        import logging

        from typing import List, Optional, Union

        import os

        from MMP import prep

        import subprocess

        from io import StringIO

        import pandas as pd

        log = logging.getLogger(__name__)



        def setup_mmp(

            csv_filename: str,

            properties_columns: Optional[List[str]] = None,

            smiles_column: str = "smiles",

            labels_column: str = "labels",

            output_filename: str = "data",

            number_of_cuts: int = 3,

            indexing_output_format: str = "mmpdb",

        ) -> int:

            """

            Function to automate from a csv file to MMP fragement and property databases

            Args:

                csv_filename (str): The file name of the input csv file containing the raw smiles strings and labels at least. This can be suplemented with properties columns.

                properties_columns (Optional[List[str]], optional): list of columns headings of properties to include or None to include zero. Defaults to None.

                smiles_column (str, optional): column heading of the smiles column. Defaults to "smiles".

                labels_column (str, optional): column heading of the label column i.e. the identifier for each molecule. Defaults to "labels".

                output_filename (str, optional): a name without and extension (i.e. no .xyz) to use to save the files to. Defaults to "data".

            Returns:

                None - saves files for future use

            """

            if output_filename.find(".") != -1:

                log.info(

                    "ouput filename given with extension. We will remove the extension and save the files using the required extensions"

                )

                output_filename = output_filename.split(".")[0].strip()

            # load the data

            raw_df = prep.load_csv_data(csv_filename)

            # Prepare the data as smiles file for fragmentation

            smi_df = prep.extract_smiles_and_labels(

                raw_df, smiles_col=smiles_column, labels_col=labels_column, valid_only=True

            )

            smi_df[[smiles_column, labels_column]].to_csv(

                f"{output_filename}.smi", sep="\t", index=False, header=None

            )

            # Make the fragment file

            _ = make_fragements(

                input_filename=f"{output_filename}.smi",

                output_filename=f"{output_filename}.fragdb",

                number_of_cuts=number_of_cuts,

            )

            # setup and create the fragment database

            _ = make_fragements_db(

                input_filename=f"{output_filename}.fragdb",

                output_filename=f"{output_filename}.{indexing_output_format}",

                output_format=f"{indexing_output_format}",

            )

            if indexing_output_format.lower().strip() == "csv":

                tmp_df = pd.read_csv(

                    f"{output_filename}.{indexing_output_format}", sep="\t", header=None

                )

                tmp_df.columns = ["smiles 1", "smiles 2", "id 1", "id 2", "v1>>v2", "constant"]

                tmp_df.to_csv(f"{output_filename}.{indexing_output_format}", index=False)

            # Prepare property data files we reload the raw data file here for safety incase there is any chance the object changed from loading before

            # TODO: determine if the data could have changed and optimize this so we can avoid reloading if the chance is low to none

            if properties_columns is not None:

                raw_df = prep.load_csv_data(csv_filename)

                properties_df = prep.extract_properties(

                    raw_df,

                    properties_columns=properties_columns,

                    smiles_col=smiles_column,

                    labels_col=labels_column,

                    valid_only=True,

                )

                properties_df.to_csv(f"properties_{output_filename}.txt", sep="\t", index=False)

                add_properties(f"properties_{output_filename}.txt", f"{output_filename}.mmpdb")

                log.info(

                    f"MMP set up completed. MMP database: {output_filename}.mmpdb with properties included"

                )

            else:

                log.info(

                    f"MMP set up completed. MMP database: {output_filename}.mmpdb without properties included"

                )



        def make_fragements(

            input_filename: str, output_filename: str, number_of_cuts: int = 3

        ) -> subprocess.CompletedProcess:

            """

            Build fragment list from smiles

            Args:

                input_filename (str): input smiles file

                output_filename (str): ouput a fragment list file

            Returns:

                subprocess.CompletedProcess: The complete process object

            Raises:

                CalledProcessError: If there is an error in running the code a raise of CalledProcessError is made with more details of issue

            """

            ret = subprocess.run(

                [

                    "mmpdb",

                    "fragment",

                    input_filename,

                    "--num-cuts",

                    str(number_of_cuts),

                    "-o",

                    output_filename,

                ],

                check=True,

                capture_output=True,

                shell=False,

            )

            if ret.returncode == 0:

                log.info(

                    f"Command completed successfully with return code {ret.returncode}. Fragment file saved as {output_filename}."

                )

                log.info(

                    f"Command output:{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            else:

                log.critical(

                    "ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully."

                )

                log.critical(

                    f"{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            return ret



        def make_fragements_db(

            input_filename: str,

            output_filename: str,

            output_format: str = "mmpdb",

        ) -> subprocess.CompletedProcess:

            """

            Build fragment database from a fragment list file (see make_fragements for the fragment list file)

            Args:

                input_filename (str): input smiles file

                output_filename (str): ouput a fragment pre database file

            Returns:

                subprocess.CompletedProcess: The complete process object

            Raises:

                CalledProcessError: If there is an error in running the code a raise of CalledProcessError is made with more details of issue

            """

            ret = subprocess.run(

                [

                    "mmpdb",

                    "index",

                    f"{input_filename}",

                    "-o",

                    f"{output_filename}",

                    "--out",

                    f"{output_format}",

                    "--symmetric",

                ],

                check=True,

                capture_output=True,

                shell=False,

            )

            if ret.returncode == 0:

                log.info(

                    f"Command completed successfully with return code {ret.returncode}. MMP database saved as {output_filename}."

                )

                log.info(

                    f"Command output:{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            else:

                log.critical(

                    "ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully."

                )

                log.critical(

                    f"{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            return ret



        def add_properties(

            input_properties_filename: str, database_name: str

        ) -> subprocess.CompletedProcess:

            """

            Add properties fields to the database which in turn enable predictions of property shifts

            Args:

                input_filename (str): input property file

                database_name (str): The name of the database to add the properties to

            Returns:

                subprocess.CompletedProcess: The complete process object

            Raises:

                CalledProcessError: If there is an error in running the code a raise of CalledProcessError is made with more details of issue

            """

            ret = subprocess.run(

                [

                    "mmpdb",

                    "loadprops",

                    database_name,

                    "--properties",

                    input_properties_filename,

                ],

                capture_output=True,

                check=True,

                shell=False,

            )

            if ret.returncode == 0:

                log.info(

                    f"Command completed successfully with return code {ret.returncode}. Properties added to MMP database {database_name} from properties file {input_properties_filename}."

                )

                log.info(

                    f"Command output:{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            else:

                log.critical(

                    "ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully."

                )

                log.critical(

                    f"{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            return ret



        def get_properties_input(

            properties_list: Union[str, List[str], None],

            property_flag: str = "--property",

            no_property_flag: str = "--no-properties",

        ) -> List[str]:

            """

            get the property choice commandline arguments

            Args:

                properties_list (Union[str, List[str]]): The properties to predict the shift in values of if any.

                The properties should be passed either as a python list or comma separated string (i.e. 'MW,Logp,LogS,Natoms').

                Alternatively pass 'none' to predict no properties or use the default 'all' to predict all known properties.

                Defaults to None which means predict all properties.

                property_flag (str, optional): The command line flag to specify a property. Defaults to "--property".

                no_property_flag (str, optional): The command line flag to not predict properties. Defaults to "--no-properties".

            Returns:

                List[str]: The commandline arguments and variables for specifitying the properties to predict if any as a python list

            """

            # determine properties to predict or none

            # First case: pass in a string 'all', 'none' (both case insenitive) or a comma separated list of properties to predict all known properties,

            # no properties or just the listed properties respectively. This interface is more for non-programatic interaction allowing strings to be

            # passed for all cases

            if isinstance(properties_list, str):

                if properties_list.lower() == "all":

                    log.info("Predicting all known properties")

                    props = list()

                elif properties_list.lower() == "none":

                    log.info("Predicting no known properties")

                    props = [no_property_flag]

                elif properties_list.find(",") != -1:

                    log.info("Predicting user selected known properties")

                    properties_list = properties_list.split(",")

                    props = []

                    for prop in properties_list:

                        props.append(property_flag)

                        props.append(prop)

            # Second case: pass in the python keyword None to make no property predictions

            elif properties_list is None:

                log.info("Predicting no known properties")

                props = [no_property_flag]

            # Third case: pass in a python iterable of properties

            else:

                log.info("Predicting user selected known properties")

                props = []

                for prop in properties_list:

                    props.append(property_flag)

                    props.append(prop)

            return props



        # TODO: this could be more generic with **kwargs test if this better or not

        def get_modifiers(

            min_variable_size: Optional[int] = None,

            max_variable_size: Optional[int] = None,

            min_constant_size: Optional[int] = None,

            max_constant_size: Optional[int] = None,

            min_radius: Optional[int] = None,

            min_pairs: Optional[int] = None,

            substructure: Optional[str] = None,

            explain: bool = False,

        ) -> List:

            """

            Prepare arguments for the transformation function as modifications of the default method

            Args:

                min_variable_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 0.

                max_variable_size (Optional[int], optional): The maximumsize in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 999.

                min_constant_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 0.

                max_constant_size (Optional[int], optional): The miaximum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 9999.

                min_radius (Optional[int], optional): The minimum environment radius around a connection site (in number of bonds) which should match for a connection to be made. Defaults to None meaning 0.

                min_pairs (Optional[int], optional): The minimum number of matched pairs (support) which needs to be reached to consider the use of a MMP rule. Defaults to None meaning 0.

                substructure (Optional[str], optional): A fixed substructure written in SMARTS which must be present to consider the transform valid. Defaults to None.

                explain (bool, optional): Provided limited SMARTS based explanations for the transformations. Defaults to False**kwargs.

            Returns:

                List: _description_

            """

            # The arguments to this function are the same names as for the commandline of MMPDB with - changed to _ for python syntax

            # In this function if any of these modifiers are set we prepare the commandline version by prepending -- and changing _ to - in the rest of the argument name

            # Must be first otherwise stores other local variables

            inps = locals().items()

            modifiers = []

            for k, v in inps:

                if v is not None and v is not False:

                    if not isinstance(k, str):

                        k = str(k)

                    arg = f"--{k.replace('_', '-')}"

                    modifiers.append(arg)

                    modifiers.append(str(v))

            return modifiers



        def generate_and_predict(

            smiles: str,

            database: str,

            properties_list: Union[str, List[str], None] = "all",

            min_variable_size: Optional[int] = None,

            max_variable_size: Optional[int] = None,

            min_constant_size: Optional[int] = None,

            max_constant_size: Optional[int] = None,

            min_radius: Optional[int] = None,

            min_pairs: Optional[int] = None,

            substructure: Optional[str] = None,

            explain: bool = False,

            **kwargs,

        ) -> pd.DataFrame:

            """

            Generate new smiles and predict the shift in the properties requested. Note the properties must exist in the database already. Please

            note that the models used to predict the shift are based on group contribution where the changes are large from the base molecules i.e.

            the known molecule these predictions can be very poor.

            Args:

                smiles (str): The smiles string base to transform and identify new possible molecules

                database (str): The MMP data base with reaction based conversion rules and properties

                properties_list (Union[str, List[str]]): The properties to predict the shift in values of if any.

                The properties should be passed either as a python list or comma separated string (i.e. 'MW,Logp,LogS,Natoms').

                Alternatively pass 'none' to predict no properties or use the default 'all' to predict all known properties.

                Defaults to None which means predict all properties.

                min_variable_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 0.

                max_variable_size (Optional[int], optional): The maximumsize in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 999.

                min_constant_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 0.

                max_constant_size (Optional[int], optional): The miaximum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 9999.

                min_radius (Optional[int], optional): The minimum environment radius around a connection site (in number of bonds) which should match for a connection to be made. Defaults to None meaning 0.

                min_pairs (Optional[int], optional): The minimum number of matched pairs (support) which needs to be reached to consider the use of a MMP rule. Defaults to None meaning 0.

                substructure (Optional[str], optional): A fixed substructure written in SMARTS which must be present to consider the transform valid. Defaults to None.

                explain (bool, optional): Provided limited SMARTS based explanations for the transformations. Defaults to False**kwargs.

            Returns:

                pd.DataFrame: The data frame of the generated molecules and predicted shifts

            """

            props = get_properties_input(properties_list=properties_list)

            modifiers = get_modifiers(

                min_variable_size,

                max_variable_size,

                min_constant_size,

                max_constant_size,

                min_radius,

                min_pairs,

                substructure,

                explain,

            )

            log.info(

                f"{' '.join(['mmpdb', 'transform', database, '--smiles', smiles, *modifiers, *props])}"

            )

            ret = subprocess.run(

                ["mmpdb", "transform", database, "--smiles", smiles, *modifiers, *props],

                check=True,

                capture_output=True,

                shell=False,

            )

            if ret.returncode == 0:

                log.info(f"Command completed successfully with return code {ret.returncode}")

                log.info(

                    f"Command output:{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            else:

                log.critical(

                    "ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully."

                )

                log.critical(

                    f"{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            transform_df = pd.read_csv(StringIO(ret.stdout.decode()), sep="\t")

            return transform_df



        # def generate()

        #     ret = subprocess.run(["mmpdb", "transform", "data.mmpdb", "--smiles", new_smiles, "--property", "prop1(madeup)"], check=True, capture_output = True, shell=False)

        #     if ret.returncode == 0:

        #         log.info(f"Command completed successfully with return code {ret.returncode}")

        #     else:

        #         log.critical("ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully.")

        # def predict()

        #     ret = subprocess.run(["mmpdb", "transform", "data.mmpdb", "--smiles", new_smiles, "--property", "prop1(madeup)"], check=True, capture_output = True, shell=False)

        #     if ret.returncode == 0:

        #         log.info(f"Command completed successfully with return code {ret.returncode}")

        #     else:

        #         log.critical("ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully.")

## Variables

```python3
log
```

## Functions


### add_properties

```python3
def add_properties(
    input_properties_filename: str,
    database_name: str
) -> subprocess.CompletedProcess
```

Add properties fields to the database which in turn enable predictions of property shifts

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| input_filename | str | input property file | None |
| database_name | str | The name of the database to add the properties to | None |

**Returns:**

| Type | Description |
|---|---|
| subprocess.CompletedProcess | The complete process object |

**Raises:**

| Type | Description |
|---|---|
| CalledProcessError | If there is an error in running the code a raise of CalledProcessError is made with more details of issue |

??? example "View Source"
        def add_properties(

            input_properties_filename: str, database_name: str

        ) -> subprocess.CompletedProcess:

            """

            Add properties fields to the database which in turn enable predictions of property shifts

            Args:

                input_filename (str): input property file

                database_name (str): The name of the database to add the properties to

            Returns:

                subprocess.CompletedProcess: The complete process object

            Raises:

                CalledProcessError: If there is an error in running the code a raise of CalledProcessError is made with more details of issue

            """

            ret = subprocess.run(

                [

                    "mmpdb",

                    "loadprops",

                    database_name,

                    "--properties",

                    input_properties_filename,

                ],

                capture_output=True,

                check=True,

                shell=False,

            )

            if ret.returncode == 0:

                log.info(

                    f"Command completed successfully with return code {ret.returncode}. Properties added to MMP database {database_name} from properties file {input_properties_filename}."

                )

                log.info(

                    f"Command output:{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            else:

                log.critical(

                    "ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully."

                )

                log.critical(

                    f"{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            return ret


### generate_and_predict

```python3
def generate_and_predict(
    smiles: str,
    database: str,
    properties_list: Union[str, List[str], NoneType] = 'all',
    min_variable_size: Optional[int] = None,
    max_variable_size: Optional[int] = None,
    min_constant_size: Optional[int] = None,
    max_constant_size: Optional[int] = None,
    min_radius: Optional[int] = None,
    min_pairs: Optional[int] = None,
    substructure: Optional[str] = None,
    explain: bool = False,
    **kwargs
) -> pandas.core.frame.DataFrame
```

Generate new smiles and predict the shift in the properties requested. Note the properties must exist in the database already. Please

note that the models used to predict the shift are based on group contribution where the changes are large from the base molecules i.e.
the known molecule these predictions can be very poor.

Args:
    smiles (str): The smiles string base to transform and identify new possible molecules
    database (str): The MMP data base with reaction based conversion rules and properties
    properties_list (Union[str, List[str]]): The properties to predict the shift in values of if any.
    The properties should be passed either as a python list or comma separated string (i.e. 'MW,Logp,LogS,Natoms').
    Alternatively pass 'none' to predict no properties or use the default 'all' to predict all known properties.
    Defaults to None which means predict all properties.
    min_variable_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 0.
    max_variable_size (Optional[int], optional): The maximumsize in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 999.
    min_constant_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 0.
    max_constant_size (Optional[int], optional): The miaximum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 9999.
    min_radius (Optional[int], optional): The minimum environment radius around a connection site (in number of bonds) which should match for a connection to be made. Defaults to None meaning 0.
    min_pairs (Optional[int], optional): The minimum number of matched pairs (support) which needs to be reached to consider the use of a MMP rule. Defaults to None meaning 0.
    substructure (Optional[str], optional): A fixed substructure written in SMARTS which must be present to consider the transform valid. Defaults to None.
    explain (bool, optional): Provided limited SMARTS based explanations for the transformations. Defaults to False**kwargs.

Returns:
    pd.DataFrame: The data frame of the generated molecules and predicted shifts

??? example "View Source"
        def generate_and_predict(

            smiles: str,

            database: str,

            properties_list: Union[str, List[str], None] = "all",

            min_variable_size: Optional[int] = None,

            max_variable_size: Optional[int] = None,

            min_constant_size: Optional[int] = None,

            max_constant_size: Optional[int] = None,

            min_radius: Optional[int] = None,

            min_pairs: Optional[int] = None,

            substructure: Optional[str] = None,

            explain: bool = False,

            **kwargs,

        ) -> pd.DataFrame:

            """

            Generate new smiles and predict the shift in the properties requested. Note the properties must exist in the database already. Please

            note that the models used to predict the shift are based on group contribution where the changes are large from the base molecules i.e.

            the known molecule these predictions can be very poor.

            Args:

                smiles (str): The smiles string base to transform and identify new possible molecules

                database (str): The MMP data base with reaction based conversion rules and properties

                properties_list (Union[str, List[str]]): The properties to predict the shift in values of if any.

                The properties should be passed either as a python list or comma separated string (i.e. 'MW,Logp,LogS,Natoms').

                Alternatively pass 'none' to predict no properties or use the default 'all' to predict all known properties.

                Defaults to None which means predict all properties.

                min_variable_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 0.

                max_variable_size (Optional[int], optional): The maximumsize in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 999.

                min_constant_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 0.

                max_constant_size (Optional[int], optional): The miaximum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 9999.

                min_radius (Optional[int], optional): The minimum environment radius around a connection site (in number of bonds) which should match for a connection to be made. Defaults to None meaning 0.

                min_pairs (Optional[int], optional): The minimum number of matched pairs (support) which needs to be reached to consider the use of a MMP rule. Defaults to None meaning 0.

                substructure (Optional[str], optional): A fixed substructure written in SMARTS which must be present to consider the transform valid. Defaults to None.

                explain (bool, optional): Provided limited SMARTS based explanations for the transformations. Defaults to False**kwargs.

            Returns:

                pd.DataFrame: The data frame of the generated molecules and predicted shifts

            """

            props = get_properties_input(properties_list=properties_list)

            modifiers = get_modifiers(

                min_variable_size,

                max_variable_size,

                min_constant_size,

                max_constant_size,

                min_radius,

                min_pairs,

                substructure,

                explain,

            )

            log.info(

                f"{' '.join(['mmpdb', 'transform', database, '--smiles', smiles, *modifiers, *props])}"

            )

            ret = subprocess.run(

                ["mmpdb", "transform", database, "--smiles", smiles, *modifiers, *props],

                check=True,

                capture_output=True,

                shell=False,

            )

            if ret.returncode == 0:

                log.info(f"Command completed successfully with return code {ret.returncode}")

                log.info(

                    f"Command output:{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            else:

                log.critical(

                    "ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully."

                )

                log.critical(

                    f"{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            transform_df = pd.read_csv(StringIO(ret.stdout.decode()), sep="\t")

            return transform_df


### get_modifiers

```python3
def get_modifiers(
    min_variable_size: Optional[int] = None,
    max_variable_size: Optional[int] = None,
    min_constant_size: Optional[int] = None,
    max_constant_size: Optional[int] = None,
    min_radius: Optional[int] = None,
    min_pairs: Optional[int] = None,
    substructure: Optional[str] = None,
    explain: bool = False
) -> List
```

Prepare arguments for the transformation function as modifications of the default method

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| min_variable_size | Optional[int] | The minimum size in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 0. | None meaning 0 |
| max_variable_size | Optional[int] | The maximumsize in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 999. | None meaning 999 |
| min_constant_size | Optional[int] | The minimum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 0. | None meaning 0 |
| max_constant_size | Optional[int] | The miaximum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 9999. | None meaning 9999 |
| min_radius | Optional[int] | The minimum environment radius around a connection site (in number of bonds) which should match for a connection to be made. Defaults to None meaning 0. | None meaning 0 |
| min_pairs | Optional[int] | The minimum number of matched pairs (support) which needs to be reached to consider the use of a MMP rule. Defaults to None meaning 0. | None meaning 0 |
| substructure | Optional[str] | A fixed substructure written in SMARTS which must be present to consider the transform valid. Defaults to None. | None |
| explain | bool | Provided limited SMARTS based explanations for the transformations. Defaults to False**kwargs. | False**kwargs |

**Returns:**

| Type | Description |
|---|---|
| List | _description_ |

??? example "View Source"
        def get_modifiers(

            min_variable_size: Optional[int] = None,

            max_variable_size: Optional[int] = None,

            min_constant_size: Optional[int] = None,

            max_constant_size: Optional[int] = None,

            min_radius: Optional[int] = None,

            min_pairs: Optional[int] = None,

            substructure: Optional[str] = None,

            explain: bool = False,

        ) -> List:

            """

            Prepare arguments for the transformation function as modifications of the default method

            Args:

                min_variable_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 0.

                max_variable_size (Optional[int], optional): The maximumsize in number of heavy atoms (non-hydrogen) to be added to a scaffold. Defaults to None meaning 999.

                min_constant_size (Optional[int], optional): The minimum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 0.

                max_constant_size (Optional[int], optional): The miaximum size in number of heavy atoms (non-hydrogen) for a fragment to be considered a scaffold. Defaults to None meaning 9999.

                min_radius (Optional[int], optional): The minimum environment radius around a connection site (in number of bonds) which should match for a connection to be made. Defaults to None meaning 0.

                min_pairs (Optional[int], optional): The minimum number of matched pairs (support) which needs to be reached to consider the use of a MMP rule. Defaults to None meaning 0.

                substructure (Optional[str], optional): A fixed substructure written in SMARTS which must be present to consider the transform valid. Defaults to None.

                explain (bool, optional): Provided limited SMARTS based explanations for the transformations. Defaults to False**kwargs.

            Returns:

                List: _description_

            """

            # The arguments to this function are the same names as for the commandline of MMPDB with - changed to _ for python syntax

            # In this function if any of these modifiers are set we prepare the commandline version by prepending -- and changing _ to - in the rest of the argument name

            # Must be first otherwise stores other local variables

            inps = locals().items()

            modifiers = []

            for k, v in inps:

                if v is not None and v is not False:

                    if not isinstance(k, str):

                        k = str(k)

                    arg = f"--{k.replace('_', '-')}"

                    modifiers.append(arg)

                    modifiers.append(str(v))

            return modifiers


### get_properties_input

```python3
def get_properties_input(
    properties_list: Union[str, List[str], NoneType],
    property_flag: str = '--property',
    no_property_flag: str = '--no-properties'
) -> List[str]
```

get the property choice commandline arguments

Args:
    properties_list (Union[str, List[str]]): The properties to predict the shift in values of if any.
    The properties should be passed either as a python list or comma separated string (i.e. 'MW,Logp,LogS,Natoms').
    Alternatively pass 'none' to predict no properties or use the default 'all' to predict all known properties.
    Defaults to None which means predict all properties.
    property_flag (str, optional): The command line flag to specify a property. Defaults to "--property".
    no_property_flag (str, optional): The command line flag to not predict properties. Defaults to "--no-properties".

Returns:
    List[str]: The commandline arguments and variables for specifitying the properties to predict if any as a python list

??? example "View Source"
        def get_properties_input(

            properties_list: Union[str, List[str], None],

            property_flag: str = "--property",

            no_property_flag: str = "--no-properties",

        ) -> List[str]:

            """

            get the property choice commandline arguments

            Args:

                properties_list (Union[str, List[str]]): The properties to predict the shift in values of if any.

                The properties should be passed either as a python list or comma separated string (i.e. 'MW,Logp,LogS,Natoms').

                Alternatively pass 'none' to predict no properties or use the default 'all' to predict all known properties.

                Defaults to None which means predict all properties.

                property_flag (str, optional): The command line flag to specify a property. Defaults to "--property".

                no_property_flag (str, optional): The command line flag to not predict properties. Defaults to "--no-properties".

            Returns:

                List[str]: The commandline arguments and variables for specifitying the properties to predict if any as a python list

            """

            # determine properties to predict or none

            # First case: pass in a string 'all', 'none' (both case insenitive) or a comma separated list of properties to predict all known properties,

            # no properties or just the listed properties respectively. This interface is more for non-programatic interaction allowing strings to be

            # passed for all cases

            if isinstance(properties_list, str):

                if properties_list.lower() == "all":

                    log.info("Predicting all known properties")

                    props = list()

                elif properties_list.lower() == "none":

                    log.info("Predicting no known properties")

                    props = [no_property_flag]

                elif properties_list.find(",") != -1:

                    log.info("Predicting user selected known properties")

                    properties_list = properties_list.split(",")

                    props = []

                    for prop in properties_list:

                        props.append(property_flag)

                        props.append(prop)

            # Second case: pass in the python keyword None to make no property predictions

            elif properties_list is None:

                log.info("Predicting no known properties")

                props = [no_property_flag]

            # Third case: pass in a python iterable of properties

            else:

                log.info("Predicting user selected known properties")

                props = []

                for prop in properties_list:

                    props.append(property_flag)

                    props.append(prop)

            return props


### make_fragements

```python3
def make_fragements(
    input_filename: str,
    output_filename: str,
    number_of_cuts: int = 3
) -> subprocess.CompletedProcess
```

Build fragment list from smiles

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| input_filename | str | input smiles file | None |
| output_filename | str | ouput a fragment list file | None |

**Returns:**

| Type | Description |
|---|---|
| subprocess.CompletedProcess | The complete process object |

**Raises:**

| Type | Description |
|---|---|
| CalledProcessError | If there is an error in running the code a raise of CalledProcessError is made with more details of issue |

??? example "View Source"
        def make_fragements(

            input_filename: str, output_filename: str, number_of_cuts: int = 3

        ) -> subprocess.CompletedProcess:

            """

            Build fragment list from smiles

            Args:

                input_filename (str): input smiles file

                output_filename (str): ouput a fragment list file

            Returns:

                subprocess.CompletedProcess: The complete process object

            Raises:

                CalledProcessError: If there is an error in running the code a raise of CalledProcessError is made with more details of issue

            """

            ret = subprocess.run(

                [

                    "mmpdb",

                    "fragment",

                    input_filename,

                    "--num-cuts",

                    str(number_of_cuts),

                    "-o",

                    output_filename,

                ],

                check=True,

                capture_output=True,

                shell=False,

            )

            if ret.returncode == 0:

                log.info(

                    f"Command completed successfully with return code {ret.returncode}. Fragment file saved as {output_filename}."

                )

                log.info(

                    f"Command output:{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            else:

                log.critical(

                    "ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully."

                )

                log.critical(

                    f"{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            return ret


### make_fragements_db

```python3
def make_fragements_db(
    input_filename: str,
    output_filename: str,
    output_format: str = 'mmpdb'
) -> subprocess.CompletedProcess
```

Build fragment database from a fragment list file (see make_fragements for the fragment list file)

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| input_filename | str | input smiles file | None |
| output_filename | str | ouput a fragment pre database file | None |

**Returns:**

| Type | Description |
|---|---|
| subprocess.CompletedProcess | The complete process object |

**Raises:**

| Type | Description |
|---|---|
| CalledProcessError | If there is an error in running the code a raise of CalledProcessError is made with more details of issue |

??? example "View Source"
        def make_fragements_db(

            input_filename: str,

            output_filename: str,

            output_format: str = "mmpdb",

        ) -> subprocess.CompletedProcess:

            """

            Build fragment database from a fragment list file (see make_fragements for the fragment list file)

            Args:

                input_filename (str): input smiles file

                output_filename (str): ouput a fragment pre database file

            Returns:

                subprocess.CompletedProcess: The complete process object

            Raises:

                CalledProcessError: If there is an error in running the code a raise of CalledProcessError is made with more details of issue

            """

            ret = subprocess.run(

                [

                    "mmpdb",

                    "index",

                    f"{input_filename}",

                    "-o",

                    f"{output_filename}",

                    "--out",

                    f"{output_format}",

                    "--symmetric",

                ],

                check=True,

                capture_output=True,

                shell=False,

            )

            if ret.returncode == 0:

                log.info(

                    f"Command completed successfully with return code {ret.returncode}. MMP database saved as {output_filename}."

                )

                log.info(

                    f"Command output:{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            else:

                log.critical(

                    "ERROR - Command returned a non-zero exit code. It has likely failed! Please check you command and input carefully."

                )

                log.critical(

                    f"{os.linesep}{ret.stdout.decode()}{os.linesep}{ret.stderr.decode()}"

                )

            return ret


### setup_mmp

```python3
def setup_mmp(
    csv_filename: str,
    properties_columns: Optional[List[str]] = None,
    smiles_column: str = 'smiles',
    labels_column: str = 'labels',
    output_filename: str = 'data',
    number_of_cuts: int = 3,
    indexing_output_format: str = 'mmpdb'
) -> int
```

Function to automate from a csv file to MMP fragement and property databases

**Parameters:**

| Name | Type | Description | Default |
|---|---|---|---|
| csv_filename | str | The file name of the input csv file containing the raw smiles strings and labels at least. This can be suplemented with properties columns. | None |
| properties_columns | Optional[List[str]] | list of columns headings of properties to include or None to include zero. Defaults to None. | None |
| smiles_column | str | column heading of the smiles column. Defaults to "smiles". | "smiles" |
| labels_column | str | column heading of the label column i.e. the identifier for each molecule. Defaults to "labels". | "labels" |
| output_filename | str | a name without and extension (i.e. no .xyz) to use to save the files to. Defaults to "data". | "data" |

**Returns:**

| Type | Description |
|---|---|
| None | None - saves files for future use |

??? example "View Source"
        def setup_mmp(

            csv_filename: str,

            properties_columns: Optional[List[str]] = None,

            smiles_column: str = "smiles",

            labels_column: str = "labels",

            output_filename: str = "data",

            number_of_cuts: int = 3,

            indexing_output_format: str = "mmpdb",

        ) -> int:

            """

            Function to automate from a csv file to MMP fragement and property databases

            Args:

                csv_filename (str): The file name of the input csv file containing the raw smiles strings and labels at least. This can be suplemented with properties columns.

                properties_columns (Optional[List[str]], optional): list of columns headings of properties to include or None to include zero. Defaults to None.

                smiles_column (str, optional): column heading of the smiles column. Defaults to "smiles".

                labels_column (str, optional): column heading of the label column i.e. the identifier for each molecule. Defaults to "labels".

                output_filename (str, optional): a name without and extension (i.e. no .xyz) to use to save the files to. Defaults to "data".

            Returns:

                None - saves files for future use

            """

            if output_filename.find(".") != -1:

                log.info(

                    "ouput filename given with extension. We will remove the extension and save the files using the required extensions"

                )

                output_filename = output_filename.split(".")[0].strip()

            # load the data

            raw_df = prep.load_csv_data(csv_filename)

            # Prepare the data as smiles file for fragmentation

            smi_df = prep.extract_smiles_and_labels(

                raw_df, smiles_col=smiles_column, labels_col=labels_column, valid_only=True

            )

            smi_df[[smiles_column, labels_column]].to_csv(

                f"{output_filename}.smi", sep="\t", index=False, header=None

            )

            # Make the fragment file

            _ = make_fragements(

                input_filename=f"{output_filename}.smi",

                output_filename=f"{output_filename}.fragdb",

                number_of_cuts=number_of_cuts,

            )

            # setup and create the fragment database

            _ = make_fragements_db(

                input_filename=f"{output_filename}.fragdb",

                output_filename=f"{output_filename}.{indexing_output_format}",

                output_format=f"{indexing_output_format}",

            )

            if indexing_output_format.lower().strip() == "csv":

                tmp_df = pd.read_csv(

                    f"{output_filename}.{indexing_output_format}", sep="\t", header=None

                )

                tmp_df.columns = ["smiles 1", "smiles 2", "id 1", "id 2", "v1>>v2", "constant"]

                tmp_df.to_csv(f"{output_filename}.{indexing_output_format}", index=False)

            # Prepare property data files we reload the raw data file here for safety incase there is any chance the object changed from loading before

            # TODO: determine if the data could have changed and optimize this so we can avoid reloading if the chance is low to none

            if properties_columns is not None:

                raw_df = prep.load_csv_data(csv_filename)

                properties_df = prep.extract_properties(

                    raw_df,

                    properties_columns=properties_columns,

                    smiles_col=smiles_column,

                    labels_col=labels_column,

                    valid_only=True,

                )

                properties_df.to_csv(f"properties_{output_filename}.txt", sep="\t", index=False)

                add_properties(f"properties_{output_filename}.txt", f"{output_filename}.mmpdb")

                log.info(

                    f"MMP set up completed. MMP database: {output_filename}.mmpdb with properties included"

                )

            else:

                log.info(

                    f"MMP set up completed. MMP database: {output_filename}.mmpdb without properties included"

                )
