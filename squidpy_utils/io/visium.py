import json

import anndata
import pandas as pd
import scipy
from skimage import io


def _read_index(fname: str):
    """Load the names from the barcode and gene name files.
    These are the output of Kallisto.
    """
    index = []
    with open(fname) as f:
        for line in f:
            stripped = line.split("\n")[0].strip()
            if len(stripped) > 0:
                index.append(stripped)
    return index


def _clean_barcode(row: pd.Series):
    """Remove the trailing -1 on the barcode names"""
    barcode = row['barcode_raw']
    barcode_clean = barcode.split('-')[0]

    return barcode_clean


def _load_unstructured_data(
        scale_factors: str,
        hires_im: str,
        lowres_im: str,
        library_id: str,
        chemistry_name: str = "Spatial 3' v1"
):
    # create unstructured data
    unstructured_data = {}

    # add the scale factors
    with open(scale_factors) as scalefactors_file:
        scalefactors_dict = json.load(scalefactors_file)
    unstructured_data.update({'scalefactors': scalefactors_dict})

    # add the images
    hires_im = io.imread(hires_im)
    lowres_im = io.imread(lowres_im)
    im_dict = {'hires': hires_im, 'lowres': lowres_im}
    unstructured_data.update({'images': im_dict})

    # add the metadata
    metadata = {'chemistry_description': chemistry_name}
    unstructured_data.update({'metadata': metadata})

    # create the final unstructured data dict
    library_id = library_id
    uns_data_dict = {'spatial': {library_id: unstructured_data}}

    return uns_data_dict


def load_visium_kallisto(
        counts_table: str,
        gene_names: str,
        barcodes: str,
        tissue_positions_list: str,
        scale_factors: str,
        hires_im: str,
        lowres_im: str,
        library_id: str,
        chemistry_name: str="Spatial 3' v1"
) -> anndata.AnnData:
    """Load a visium dataset that was preprocessed with Kallisto

    Parameters
    ----------
    counts_table : str
        The path to the counts table output from Kallisto.
        This file usually has the extension ".mtx".
    gene_names : str
        The path to the file containing the gene names for the
        Kallisto counts table. This file usually ends with "genes.txt".
    barcodes : str
        The path to the file containing the spot barcodes for the
        Kallisto counts table. This file usually ends with "barcodes.txt".
    tissue_positions_list : str
        The path to the file containing the coordinates of each barcode
        that is output from the 10X space ranger pipeline.
        This file is usually called: "tissue_positions_list.csv".
    scale_factors : str
        The path to the file output by the 10X space ranger pipeline
        containing the scale factors that map the hires and
        lowres image to the original image.
        This file is usually called: "scalefactors_json.json"
    hires_im : str
        The path to the hires image that is output from the
        10X space ranger pipeline.
        This file is usually called: tissue_hires_image.png
    lowres_im : str
        The path to the lowres image that is output from the
        10X space ranger pipeline.
        This file is usually called: tissue_lowres_image.png
    library_id  : str
        The unique identifier for the library that was sequenced.
    chemistry_name : str
        The name of the chemistry used to create the library.
        The default value is: "Spatial 3' v1"

    Returns
    -------
    adata : anndata.AnnData
        The AnnData object containing the visium results
    """
    # load the gene and barcode titles
    genes = _read_index(gene_names)
    barcodes = _read_index(barcodes)
    ordered_spot_data = pd.DataFrame({'barcode': barcodes})

    # load the count table
    bg_matrix = scipy.io.mmread(counts_table).toarray()

    # load the data about the spots (for AnnData.obs)
    spot_coords = pd.read_csv(
        tissue_positions_list,
        header=None,
        names=[
            'barcode_raw',
            'in_tissue',
            'array_col',
            'array_row',
            'im_x',
            'im_y'
        ]
    )

    # the barcode names in the tissue_position_list.csv usually
    # have a trailing -1, so we have to remove it
    spot_coords['barcode'] = spot_coords.apply(_clean_barcode, axis=1)

    # get the spots that are in the table
    spots_in_data = spot_coords.loc[
        spot_coords['barcode'].isin(ordered_spot_data['barcode'])
    ]
    ordered_spot_data = pd.merge(ordered_spot_data, spots_in_data, on='barcode')

    # AnnData makes it a string index
    ordered_spot_data.index = ordered_spot_data.index.map(str)

    unstructured_data = _load_unstructured_data(
        scale_factors=scale_factors,
        hires_im=hires_im,
        lowres_im=lowres_im,
        library_id=library_id,
        chemistry_name=chemistry_name
    )

    # construct the final AnnData object
    adata = anndata.AnnData(
        X=bg_matrix,
        obs=ordered_spot_data,
        uns=unstructured_data,
    )
    adata.obsm['spatial'] = ordered_spot_data[['im_y', 'im_x']].values
    adata.var_names = genes

    return adata
