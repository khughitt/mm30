"""
create_combined_sample_metadata.py
"""
import re
import glob
import pandas as pd

geo_mdata = glob.glob(config['sample_metadata']['geo'])
mmrf_mdata = config['sample_metadata']['mmrf'] 

# combine geo/mmrf metadata into a single dataframe with columns for sample id,
# experiment, platform, and platform_type
mdat = pd.read_feather(mmrf_mdata)[["public_id", "disease_stage"]]

mdat.rename(columns={'public_id': 'sample_id'}, inplace=True)

# mmrf
mdat['experiment'] = ['MMRF'] * mdat.shape[0]
mdat['platform_id'] = ['GPL16791'] * mdat.shape[0] # HiSeq 2500
mdat['platform_type'] = ['RNA-Seq'] * mdat.shape[0]

# microarray platforms included in MM30
mm30_microarray_platforms = ["GPL96", "GPL97", "GPL570", "GPL10558", "GPL5175", "GPL6244", "GPL25401", "GPL4819"]

# geo
for infile in geo_mdata:
    # determine GEO identifier
    geo_id = re.findall('GSE[0-9]+', infile)[0]

    # load geo metadata
    geo_mdat = pd.read_feather(infile)

    # rename sample id column
    geo_mdat.rename(columns={geo_mdat.columns[0]: 'sample_id'}, inplace=True)

    # add experiment column
    geo_mdat['experiment'] = [geo_id] * geo_mdat.shape[0]

    # determine platform type
    if geo_mdat.platform_id.iloc[0] in mm30_microarray_platforms:
        geo_mdat['platform_type'] = ['Microarray'] * geo_mdat.shape[0]
    else:
        geo_mdat['platform_type'] = ['RNA-Seq'] * geo_mdat.shape[0]

    # match column order
    try:
        geo_mdat = geo_mdat[mdat.columns]
    except:
        print(infile)
        breakpoint()

    # append to combined dataframe
    mdat = pd.concat([mdat, geo_mdat])

# write combined metadata to disk
mdat.reset_index(drop=True).to_feather(snakemake.output[0])
