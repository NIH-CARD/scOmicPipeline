from anndata import AnnData
import pickle
import argparse

"""
Allows multiple arguments per flag to be read
"""
class MultiLineArgAndFileParser(argparse.ArgumentParser): 
    def convert_arg_line_to_args(self, arg_line): 
        return arg_line.split() 

"""
Loads a pickle if no anndata object is provided
"""
def choose_adata_src(input):
    if isinstance(input, AnnData):
        adata = input
    else:
        with open(input, 'rb') as file:
            adata = pickle.load(file)

    return adata

"""
Create an input name based on the provided name 
and the project name
"""
def get_input_name(project, inp):
    return project + '_' + inp


"""
Get an output name based on the provided name
and the project name.
"""
def get_output_name(project, out):
    return project + '_' + out

"""
Dump to pickle file
"""
def dump_to_pickle(out, adata):
    # Dump adata out
    with open(out, 'wb') as file:
        pickle.dump(adata, file)


"""
Returns saving names more easily
series = name to append after project
returns full name
"""
def get_save_name(folder="", prefix="", series=None, project="", figure_type=""):
    if series:
        series = "_" + series + "_"
    else:
        series = "_"
    save_name = prefix + project + series + "." + figure_type
    
    return str(save_name)

