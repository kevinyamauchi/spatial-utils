import numpy as np
import pandas as pd

from spatial_utils.io.visium import _clean_barcode


def test_clean_barcode():
    barcodes = [str(i) for i in range(10)]
    barcodes_raw = [b + '-1' for b in barcodes]

    df = pd.DataFrame({'barcode_raw': barcodes_raw})
    df['barcode'] = df.apply(_clean_barcode, axis=1)

    np.testing.assert_array_equal(df['barcode'], barcodes)


def test_already_cleaned_barcode():
    barcodes = [str(i) for i in range(10)]

    df = pd.DataFrame({'barcode_raw': barcodes})
    df['barcode'] = df.apply(_clean_barcode, axis=1)

    np.testing.assert_array_equal(df['barcode'], barcodes)
