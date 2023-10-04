import importlib.resources as res
import pandas as pd

package = res.files('cefcon')
TFs_human = pd.read_csv(package / 'resources/hs_hgnc_tfs_lambert2018.txt', names=['gene_symbol'])
TFs_mouse = pd.read_csv(package / 'resources/mm_mgi_tfs.txt', names=['gene_symbol'])
TFs_human_animaltfdb = pd.read_csv(package / 'resources/hs_hgnc_tfs_animaltfdb4.txt', names=['gene_symbol'])
TFs_mouse_animaltfdb = pd.read_csv(package / 'resources/mm_mgi_tfs_animaltfdb4.txt', names=['gene_symbol'])

__all__ = ['TFs_human', 'TFs_mouse', 'TFs_human_animaltfdb', 'TFs_mouse_animaltfdb']
