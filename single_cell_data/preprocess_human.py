import os
import argparse
import scanpy as sc
import mygene

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", required=True, type=str, help="string, path which contains .mtx file")
    parser.add_argument("-i", "--folder", required=True, type=str, help="string, path to the folder to save the data")
    parser.add_argument("-k", "--topk", default=2000, type=int, help="integer, optional, number of most variable genes,"
                                                                     " default 2000")
    parser.add_argument("-d", "--dataset", required=True, type=str, help="choose which datasets to use")
    return parser.parse_args()


def main(args):
    if args.dataset == 'altas_lung':    # args.path : human/altas_lung/
        data = sc.read_10x_mtx(args.path, var_names="gene_symbols", cache=False)
        sc.pp.filter_cells(data, min_genes=500)
        sc.pp.filter_genes(data, min_cells=5)
        # sc.pp.highly_variable_genes(data,min_disp=0.01,min_mean=0,max_mean=1000,flavor='seurat',subset=True,inplace=True,batch_key=None,check_values=True)
        sc.pp.highly_variable_genes(data, min_disp=0.01, min_mean=0, max_mean=1000, subset=True, inplace=True,
                                    batch_key=None, check_values=True, flavor='seurat_v3', n_top_genes=args.topk)
        # We then normalize each cell to 1e4 total read counts.
        sc.pp.normalize_total(data, target_sum=1e4)
        # The values are then log transformed and scaled to zero mean and unit variance.
        sc.pp.log1p(data)
        sc.pp.scale(data, max_value=10)

    elif args.dataset == 'pbmx28k':   # args.path : human/pbmc28k/fresh_68k_pbmc_donor_a_raw_gene_bc_matrices/hg19/
        data = sc.read_10x_mtx(args.path, var_names="gene_symbols", cache=False)
        sc.pp.filter_cells(data, min_genes=200)
        sc.pp.filter_genes(data, min_cells=3)
        # sc.pp.highly_variable_genes(data, min_disp=0.01, min_mean=0, max_mean=1000, subset=True, inplace=True,
        #                             batch_key=None, check_values=True, flavor='seurat_v3', n_top_genes=args.topk)
        sc.pp.highly_variable_genes(data, flavor='seurat_v3', n_top_genes=args.topk)
        sc.pp.normalize_total(data, target_sum=1e4)
        sc.pp.log1p(data)
        sc.pp.scale(data, max_value=10)


    elif args.dataset == 'hubmap':
        data = sc.read(args.path, dtype='float64')

        # For each dataset, we first filter out cells expressing fewer than 200 genes and genes
        # expressed in fewer than 3 cells.
        sc.pp.filter_cells(data, min_genes=200)
        sc.pp.filter_genes(data, min_cells=3)
        sc.pp.highly_variable_genes(data, flavor='seurat_v3', n_top_genes=args.topk)
        # We then normalize each cell to 1e4 total read counts.
        sc.pp.normalize_total(data, target_sum=1e4)
        # The values are then log transformed and scaled to zero mean and unit variance.
        sc.pp.log1p(data)
        # We clip the values such that the maximum is ten.
        sc.pp.scale(data, max_value=10)

        mg = mygene.MyGeneInfo()
        new_idx = []
        for each in data.var_names:
            new_idx.append(each.split('.')[0])
        new_index = mg.querymany(new_idx, scopes='ensembl.gene', fields='symbol', species='human')
        new_index_map = {}
        for item in new_index:
            query = item['query']
            new_name = item['symbol'] if 'symbol' in item else item['query']
            new_index_map[query] = new_name
        data.var.rename(index=new_index_map)

    # data = sc.read_10x_mtx(args.path, var_names="gene_symbols", cache=False)
    # sc.pp.filter_cells(data, min_genes=500)
    # sc.pp.filter_genes(data, min_cells=5)
    # # We then normalize each cell to 1e4 total read counts.
    # sc.pp.normalize_total(data, target_sum=1e4)
    # # The values are then log transformed and scaled to zero mean and unit variance.
    # sc.pp.log1p(data)
    # data.write_h5ad(os.path.join(args.folder, "human_lung_processed.h5ad"))
    data.write_h5ad(os.path.join(args.folder, args.dataset+"_processed.h5ad"))



if __name__ == "__main__":
    parsed_args = parse_args()
    main(parsed_args)
