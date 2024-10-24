from typing import Optional

import os
from pathlib import Path

import pandas as pd
import numpy as np
from itertools import permutations, product
import scanpy as sc
import requests
import zipfile
from tqdm.auto import tqdm


def _download_from_url(file_url: str, save_path: Path):
    """
    Downloads a file from a given URL to a specified save path.

    Parameters:
        file_url (str): The URL from which the file needs to be downloaded.
        save_path (Path): The path where the downloaded file needs to be saved.
    """
    try:
        response = requests.get(file_url, stream=True)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f'Download error: {e}')
        return

    total_size = int(response.headers.get('content-length', 0))
    block_size = 8192
    progress_bar = tqdm(total=total_size, unit='B', unit_scale=True)

    download_file = save_path.parent / file_url.split('/')[-1]
    with open(download_file, 'wb') as file:
        for chunk in response.iter_content(chunk_size=block_size):
            if chunk:
                file.write(chunk)
                progress_bar.update(len(chunk))

    progress_bar.close()

    if str(download_file).endswith('.zip'):
        with zipfile.ZipFile(download_file, 'r') as zip_file:
            zip_file.extractall(os.path.dirname(save_path))
            zip_file.close()
        os.remove(download_file)

    print(f'Ths data has been downloaded to `{save_path}`.')


def load_human_prior_interaction_network(dataset: str = 'nichenet',
                                         only_directed: bool = False,
                                         force_download: bool = False):
    """
    Load and process a human prior gene interaction network dataset.

    Parameters:
        dataset (str): The name of the dataset to load. Default is 'nichenet'.
            Available options: {'nichenet', 'pathwaycommons', 'inbiomap', 'harmonizome', 'omnipath_interactions'}.
        only_directed (bool): Whether to only include directed edges in the network. Default is False.
        force_download (bool): Whether to force the download of the dataset, even if it already exists locally. Default is False.

    Returns:
        (pd.DataFrame): A DataFrame containing the processed prior gene interaction network.
            The DataFrame has two columns: 'from' and 'to', representing the source and target genes of the edges.
    """

    # The URL for every dataset. These datasets are stored at zenodo (https://doi.org/10.5281/zenodo.7564872).
    urls = {
        'nichenet': 'https://zenodo.org/record/8013900/files/NicheNet_human.zip',
        'pathwaycommons': 'https://zenodo.org/record/8013900/files/PathwayCommons12.All.hgnc.zip',
        'inbiomap': 'https://zenodo.org/record/8013900/files/InBioMap.zip',
        'harmonizome': 'https://zenodo.org/record/8013900/files/Harmonizome_nichenet.zip',
        'omnipath_interactions': 'https://zenodo.org/record/8013900/files/Omnipath_interaction.zip',
    }
    filenames = {
        'nichenet': 'NicheNet_human.csv',
        'pathwaycommons': 'PathwayCommons12.All.hgnc.sif',
        'inbiomap': 'InBioMap.csv',
        'harmonizome': 'Harmonizome_nichenet.csv',
        'omnipath_interactions': 'Omnipath_interaction.csv',
    }

    # Download if the file does not exist locally or 'force_download' is True
    data_path = Path('./data_cache') / filenames[dataset]
    if force_download or not data_path.exists():
        data_path.parent.mkdir(parents=True, exist_ok=True)
        _download_from_url(urls[dataset], data_path)

    # Process the dataset based on its name
    if dataset == 'nichenet':  # 5,583,023
        prior_net = pd.read_csv(data_path, index_col=None, header=0)

    elif dataset == 'pathwaycommons':  # 1,200,159
        prior_net = pd.read_csv(data_path, sep='\t',
                                names=['from', 'type', 'to'])
        undirected_type = ['interacts-with', 'in-complex-with']
        type_mapper = dict({v: 1 for v in list(set(prior_net.type.unique()) - set(undirected_type))},
                           **{v: 0 for v in undirected_type})
        prior_net['is_directed'] = prior_net['type'].map(type_mapper)
        prior_net = prior_net.loc[~prior_net['from'].str.startswith('CHEBI'), :]
        prior_net = prior_net.loc[~prior_net['to'].str.startswith('CHEBI'), :]
        del prior_net['type']

    elif dataset == 'inbiomap':  # 625,641
        prior_net = pd.read_csv(data_path)
        prior_net.rename(columns={'genesymbol_a': 'from', 'genesymbol_b': 'to'}, inplace=True)
        prior_net = prior_net.loc[:, ['from', 'to']]

    elif dataset == 'harmonizome':  # 3,418,949
        prior_net = pd.read_csv(data_path, sep='\t')
        prior_net = prior_net.loc[:, ['from', 'to']]

    elif dataset == 'omnipath_interactions':  # 525,430
        prior_net = pd.read_csv(data_path)
        prior_net.rename(columns={'source_genesymbol': 'from', 'target_genesymbol': 'to'}, inplace=True)
        complex_idx = prior_net['source'].str.startswith('COMPLEX') | \
                      prior_net['target'].str.startswith('COMPLEX')
        prior_net_genes_only = prior_net.loc[~complex_idx, ['from', 'to', 'is_directed']]
        prior_net_genes_complex = prior_net.loc[complex_idx, ['from', 'to', 'is_directed']]
        # process complex items
        temp_edge = []
        temp_edge_type = []
        for source, target in zip(prior_net_genes_complex['from'], prior_net_genes_complex['to']):
            source = source.split('_')
            target = target.split('_')
            # inside the complex
            intra = list(permutations(source, r=2)) + list(permutations(target, r=2))
            temp_edge += intra
            temp_edge_type += [0] * len(intra)
            # between two complexes
            inter = list(product(source, target))
            temp_edge += inter
            temp_edge_type += [1] * len(inter)
        prior_net_complex = pd.DataFrame(temp_edge, columns=['from', 'to'])
        prior_net_complex['is_directed'] = temp_edge_type
        prior_net = pd.concat([prior_net_genes_only, prior_net_complex], axis=0)
        prior_net = prior_net.dropna()
        prior_net = prior_net.drop_duplicates()

    else:
        print(f"Value error. {dataset} is not available.")
        print("Available option: {'nichenet', 'pathwaycommons', 'inbiomap', 'harmonizome', 'omnipath_interactions'}")

    if ('is_directed' in prior_net) and only_directed:
        prior_net = prior_net[prior_net['is_directed'] == 1]

    prior_net = prior_net[['from', 'to']].drop_duplicates().astype(str)

    print(f"Load the prior gene interaction network: {dataset}. "
          f"#Genes: {len(np.unique(prior_net.iloc[:, [0, 1]]))}, "
          f"#Edges: {len(prior_net)}")

    return prior_net


def convert_human_to_mouse_network(net: pd.DataFrame, save: bool = False):
    """
    Converts the gene symbols of a human prior gene interaction network to mouse gene symbols.

    Parameters:
        net (pd.DataFrame): A DataFrame representing the human prior gene interaction network.
                            It should have two columns named 'from' and 'to', representing the interacting genes.
        save (bool): A flag indicating whether to save the converted network. Default is False.

    Returns:
        (pd.DataFrame): A DataFrame representing the converted mouse prior gene interaction network.
    """
    import biomart

    print('Convert genes of the prior gene interaction network to mouse gene symbols.')
    print('Linking the Ensembl server...')
    with tqdm(total=10, desc='Processing', miniters=1) as outer_bar:
        outer_bar.update()

        # Set up connection to server
        for name in ['ensembldb', 'asia', 'useast', 'martdb']:
            try:
                server = biomart.BiomartServer(
                    f'http://{name}.ensembl.org/biomart/')
                print(f'Server \'http://{name}.ensembl.org/biomart/\' is OK')
                break
            except Exception as e:
                print(f'[Error] Server not available: http://{name}.ensembl.org/biomart//martservice')

        human_dataset = server.datasets['hsapiens_gene_ensembl']
        outer_bar.update()
        mouse_dataset = server.datasets['mmusculus_gene_ensembl']
        outer_bar.update()

        human_attributes = ['ensembl_gene_id', 'hgnc_symbol']
        mouse_attributes = ['ensembl_gene_id', 'mgi_symbol']  # 'external_gene_name'
        to_homolog_attribute = 'mmusculus_homolog_ensembl_gene'

        # Map gene symbol to ensembl ID of query species
        query = human_dataset.search({'attributes': human_attributes})
        query = query.raw.data.decode('ascii').split('\n')[:-1]
        query = pd.DataFrame([d.split('\t') for d in query], columns=['human_ensembl_id', 'hgnc_symbol'])
        outer_bar.update(2)

        # Map ensembl IDs between two species
        from2to = human_dataset.search({'attributes': ['ensembl_gene_id', to_homolog_attribute]})
        from2to = from2to.raw.data.decode('ascii').split('\n')[:-1]
        from2to = pd.DataFrame([d.split('\t') for d in from2to], columns=['human_ensembl_id', 'mouse_ensembl_id'])
        from2to = from2to.merge(query, how='outer', on='human_ensembl_id')
        outer_bar.update()

        # Map ensembl ID to gene symbol of target species
        target = mouse_dataset.search({'attributes': mouse_attributes})
        target = target.raw.data.decode('ascii').split('\n')[:-1]
        target = pd.DataFrame([d.split('\t') for d in target], columns=['mouse_ensembl_id', 'mgi_symbol'])
        target = target.merge(from2to, how='outer', on='mouse_ensembl_id')
        outer_bar.update()

        # Gene mapper of the network
        query_genes = np.unique(net.loc[:, ['from', 'to']].astype(str))
        mapper = target.loc[target['hgnc_symbol'].isin(query_genes), ['hgnc_symbol', 'mgi_symbol']].copy()
        mapper = mapper.dropna()
        mapper = mapper.drop_duplicates()
        mapper = mapper.loc[mapper['mgi_symbol'] != '', :]

        # Process ambiguous (1-to-many) and unambiguous (1-to-1 and many-to-1) genes separately
        human_gene_value_counts = mapper.loc[:, 'hgnc_symbol'].value_counts()
        unambiguous_genes = human_gene_value_counts[human_gene_value_counts == 1].index.tolist()
        ambiguous_genes = human_gene_value_counts[human_gene_value_counts > 1].index.tolist()
        outer_bar.update()

        # Directly convert the interactions with unambiguous genes
        net_una = net.loc[net['from'].isin(unambiguous_genes) & net['to'].isin(unambiguous_genes), ['from', 'to']]
        converted_network_unambiguous = pd.DataFrame()
        mapper_una_dict = mapper.loc[mapper['hgnc_symbol'].isin(unambiguous_genes), :].set_index(['hgnc_symbol'])[
            'mgi_symbol'].to_dict()
        converted_network_unambiguous['from'] = net_una['from'].map(mapper_una_dict)
        converted_network_unambiguous['to'] = net_una['to'].map(mapper_una_dict)
        outer_bar.update()

        # Process interactions where one gene is ambiguous and another is unambiguous
        net_a = net.loc[(net['from'].isin(ambiguous_genes) & net['to'].isin(unambiguous_genes)) |
                        (net['to'].isin(ambiguous_genes) & net['from'].isin(unambiguous_genes)), ['from', 'to']]
        mapper_a = mapper.loc[mapper['hgnc_symbol'].isin(ambiguous_genes), :]
        temp_edge = []
        with tqdm(total=len(net_a),
                  desc='Converting ambiguous gene symbols',
                  leave=False,
                  miniters=1,
                  ) as pbar:
            for source, target in zip(net_a['from'], net_a['to']):
                if source in unambiguous_genes:
                    source_convert = [mapper_una_dict[source]]
                    target_convert = mapper_a[mapper_a['hgnc_symbol'] == target]['mgi_symbol'].tolist()
                else:
                    source_convert = mapper_a[mapper_a['hgnc_symbol'] == source]['mgi_symbol'].tolist()
                    target_convert = [mapper_una_dict[target]]
                temp_edge += list(product(source_convert, target_convert))
                pbar.update(1)
        converted_network_ambiguous = pd.DataFrame(temp_edge, columns=['from', 'to'])
        outer_bar.update()

    # Combine the converted network
    prior_net_converted = pd.concat([converted_network_unambiguous, converted_network_ambiguous], axis=0)
    prior_net_converted = prior_net_converted.drop_duplicates()

    print(f"The converted prior gene interaction network: "
          f"#Genes: {len(np.unique(prior_net_converted.iloc[:, [0, 1]]))}, "
          f"#Edges: {len(prior_net_converted)}")

    if save:
        file_name = Path('./data_cache') / 'prior_net_for_mouse.csv'
        prior_net_converted.to_csv(file_name)
        print(f'Ths converted network has been save to `{file_name}`.')

    return prior_net_converted


def mouse_hsc_nestorowa16(fpath: Optional[str] = './data_cache/mouse_hsc_nestorowa16_v0.h5ad', version: Optional[str] = 'v0'):
    if version=='v0':
        fpath = './data_cache/mouse_hsc_nestorowa16_v0.h5ad'
        url = 'https://zenodo.org/record/8013900/files/mouse_hsc_nestorowa16_v0.h5ad'
        print('Load mouse_hsc_nestorowa16_v0.h5ad')
    elif version=='v1':
        fpath = './data_cache/mouse_hsc_nestorowa16_v1.h5ad'
        url = 'https://zenodo.org/record/8013900/files/mouse_hsc_nestorowa16_v1.h5ad'
        print('Load mouse_hsc_nestorowa16_v1.h5ad')
    else:
        raise ValueError("Value error. version should be 'v0' or 'v1'.")

    adata = sc.read(fpath, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


def mouse_hsc_paul15(fpath: Optional[str] = './data_cache/mouse_hsc_paul15.h5ad'):
    url = 'https://zenodo.org/record/7564872/files/mouse_hsc_paul15.h5ad'
    adata = sc.read(fpath, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata