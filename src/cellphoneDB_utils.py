import os
import sys
import json
import random
import logging
import networkx as nx
import scipy.stats
import pandas as pd
import numpy as np
import seaborn as sns
from itertools import *
import matplotlib.pyplot as plt
from multiprocessing import Pool, Manager
from scipy.cluster.hierarchy import ward
from sklearn.neighbors import KernelDensity

logger = logging.getLogger("cellphoneDB_utils")
logger.setLevel(logging.INFO)
handler=logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s:%(step)s - %(levelname)s - %(message)s"))
logger.addHandler(handler)

def process_metadata(data_input):
    """Transform cell type long dataframe to wider table with original cell index
    Parameters:
    ------------
    data_input: dataframe
        Includes default column "cell_type" and index "cell id"
    col_cell_type: string
        Column name of cell type,  cell type names must be syntactically valid

    Returns:
    data_out: dataframe
        dummy data frame with cell types on columns
    """

    original_index = data_input['id_new']
    data_input["value"] = 1
    data_out = data_input[["id_new", "cell_type", "value"]].pivot_table(index="id_new", columns="cell_type", values="value",fill_value=0).reindex(original_index)

    return data_out

def sp_grid_kern_bin(data, coord, min_num=15, h=20, n=100j, tot_num=True):
    """ construct the map of grid from data and return Kenerl density on grid

    Parameters
    ---------
    data: dataframe
        Cell type dummy table
    coord: dataframe
        Coordinates data with X and Y columns
    min_num: int
        Minimum cells count for each cell type
    h: int
        Bandwidths for x and y directions, more details in KernelDensity function in sklearn.neighbors package
    n: int + j
        Number of grid points in each direction
    tot_num: bool to decide whether to normalize

    Returns:
    ---------
    data_out: dataframe with X1,X2 coords and kernel density
        data_out with cell types columns representing kernel density in each spot.

    """

    # for cells with less than min_num, randomly add 1
    data_ = data
    if min_num < 0: min_num = 1
    col_min = data.sum() < min_num
    col_min_index = col_min[col_min].index
    if len(col_min_index) != 0:
        for i in range(0, len(col_min_index)):
            random.seed(i)
            cell_type = col_min_index[i]
            ct_dummy_array = data_.loc[:, cell_type]
            random_cell = random.sample(population=list(np.where(ct_dummy_array == 0)[0]),
                                        k=min_num - ct_dummy_array.sum())
            data_.loc[:, cell_type][random_cell] = 1
            print(col_min_index[i] + " cell randomly add to " + str(data_.loc[:, col_min_index[i]].sum()) + " cells")

    # coordinates for X1 and X2
    coord = pd.DataFrame(coord)
    coord.columns = ['X' + str(i) for i in range(0, len(coord.columns))]
    coord.index = list(data_.index)
    data_merge = pd.concat([coord, data_], axis=1)

    kde2d = KernelDensity(bandwidth=h, kernel="gaussian")
    # build 2D grid as sample points
    xx, yy = np.mgrid[data_merge['X0'].min():data_merge['X0'].max():n,
             data_merge['X1'].min():data_merge['X1'].max():n]
    xy_sample = pd.DataFrame(np.vstack([xx.ravel(), yy.ravel()]).T, columns=["X", "Y"])
    data_out = pd.DataFrame(xy_sample)
    data_out.columns = ["X0", "X1"]

    # print("estimate gaussian kernel 2D density for each cell types...")
    for i in data_.columns:
        xy_train = data_merge[["X0", "X1"]][data_merge[i] == 1]
        kde2d.fit(xy_train)
        # score_samples() returns the log-likelihood of the samples
        z = np.exp(kde2d.score_samples(xy_sample))

        # plt.figure(figsize=(6,8))
        # plt.pcolormesh(xx, yy, np.reshape(z, xx.shape), cmap="YlGnBu")
        # plt.scatter(xy_train['X0'], xy_train['X1'], s=2, facecolor="white")
        # plt.savefig("Cell type gaussian kernel density/" +i+".pdf")
        z = pd.DataFrame(z, columns=[i])
        data_out = pd.concat([data_out, z], axis=1)

    # data_out with cell types columns representing kernel density in each spot. (normalization)
    if tot_num:
        data_out = data_out.drop(["X0", "X1"], axis=1)
        data_out = data_out / data_out.sum()

    return data_out

def KL_JS_Divergence(X, eps=1e-20, diver="KL"):
    """ calculate Kullback-Leibler or Jensen-Shannon diversity

    Parameters
    ---------
    X: dataframe
        Density matrix
    eps: double flout
        small value added
    diver: str
        "KL" or "JS" representing Kullback-Leibler or Jensen-Shannon Divergence

    Returns
    ---------
    KL_D: dataframe
        KL-divergence matrix

    """

    X_ = X.fillna(eps)
    X_[X_ < eps] = eps
    X_ = X_ / X_.sum()
    n_type = len(X_.columns)
    diver_matrix = pd.DataFrame(np.zeros((n_type, n_type)))
    diver_matrix.index = X_.columns
    diver_matrix.columns = X_.columns
    # print("calculate cell types pairs " + diver + " divergence...")
    if diver == "KL":
        for i in combinations(X_.columns, 2):
            KL = scipy.stats.entropy(X_[i[0]], X_[i[1]])
            diver_matrix.loc[i[0], i[1]] = KL
            diver_matrix.loc[i[1], i[0]] = KL
    else:
        for i in combinations(X_.columns, 2):
            M = (X_[i[0]] + X_[i[1]]) / 2
            JS = 0.5 * scipy.stats.entropy(X_[i[0]], M, base=2) + 0.5 * scipy.stats.entropy(X_[i[1]], M, base=2)
            diver_matrix.loc[i[0], i[1]] = JS
            diver_matrix.loc[i[1], i[0]] = JS

    return diver_matrix

def boot(res_dict, i, boot_n, data_boot, coord_boot, min_num, h, tot_num, diver):
    print('---- Bootstrap ' + str(i + 1) + " time Begin ----")
    k2d_boot = sp_grid_kern_bin(data=data_boot, coord=coord_boot, min_num=min_num, h=h, tot_num=tot_num)
    dis_boot = KL_JS_Divergence(k2d_boot, eps=1e-20, diver=diver)

    # create a graph from the adjacency matrix
    graph_boot = nx.from_pandas_adjacency(-np.log2(dis_boot))
    # MST
    graph_mst_boot = nx.maximum_spanning_tree(graph_boot)
    mst_boot = nx.to_pandas_adjacency(graph_mst_boot)
    # cumulate bootstrap
	## 每轮bootstrap所有的结果存储
    res_dict["dis_boot_array"][:, :, i] = dis_boot
    ## 每轮bootstrap所有的结果的平均值存储
    res_dict["dis_cons"] = res_dict["dis_cons"] + dis_boot / boot_n
    ## 每轮bootstrap所有的最小生成树结果的平均结果
    res_dict["mst_cons"] = res_dict["mst_cons"] + mst_boot / boot_n
    print('---- Bootstrap ' + str(i + 1) + " time End ----")


def KL_JS_boot_mst(dummy_df, coord_df, min_num=15, boot_n=20, prop=0.8, h=20, tot_num=True, diver="JS",n_threads=30):
    """  calculate KL or JS divergence and use MST to generate a tree structure by Bootstrap
         to obtain consensus cell types colocalization dissimilarity matrix

    Parameters
    ---------
    dummy_df: dataframe
        dummy data frame with cell types on columns
    coord_df: dataframe
        Coordinates data with X and Y columns
    min_num: int
        Minimum cells count for each cell type
    boot_n: int
        Number of bootstraping iteration
    prop: 0-1
        Subsample preportion
    diver: String
        use KL or JS divergence

    other see in sp_grid_kern_bin function

    Returns
    ---------
    dis_cons: dataframe
        Bootstrap KL/JS divergence
    mst_cons: dataframe
        Bootstrap MST matrix
    dis_boot_array: array
        Each Bootstrap result

    """
    sys.stdout.write('### Begin Predicting Spatially Co-localized Cell Types by Bootstrap ###' + "\n")
    sys.stdout.write('### Methods: ' + diver + ', MST ###' + "\n")
    n_smp = len(dummy_df)
    n_type = len(dummy_df.columns)
    # Process Pool for maximum 3 processes
    pool = Pool(processes=n_threads)
    # shared data
    manager = Manager()
    res_dict = manager.dict()
    res_dict["dis_boot_array"] = np.zeros((n_type, n_type, boot_n))
    # dis_cons
    res_dict["dis_cons"] = pd.DataFrame(np.zeros((n_type, n_type)), columns=list(dummy_df.columns),
                                        index=list(dummy_df.columns))
    # mst_cons
    res_dict["mst_cons"] = pd.DataFrame(np.zeros((n_type, n_type)), columns=list(dummy_df.columns),
                                        index=list(dummy_df.columns))
    for i in range(boot_n):
        random.seed(i)
        # Bootstrap
				# 从n_smp长度的例表中随机抽取80%的样本
        idx = random.sample(range(n_smp), round(n_smp * prop))
        data_boot = dummy_df.iloc[idx, :]
        coord_boot = coord_df.iloc[idx, :]
        pool.apply_async(boot, (res_dict, i, boot_n, data_boot, coord_boot, min_num, h, tot_num, diver))
    # 扔了 1000个进程进进程池后，关闭进程池，不允许新的进程加入
    pool.close()
    # 运行进程池中的进程
    pool.join()

    return res_dict["dis_boot_array"], res_dict["dis_cons"], res_dict["mst_cons"]

## 通过最小生成树的结果卡阈值，然后定义微环境，保存为cellphone可以读入的格式文件
def network_microenv(df_adjacency, out_path, to_cpdb=True, cutoff=0.5):
## df_adjacency:每轮bootstrap所有的最小生成树结果的平均结果(mst_cons)
    """ obtain cell types microenvironment

        Parameters
        ---------
        df_adjacency: dataframe / matrix
            divergence matrix
        out_path: string
            picture name
        cutoff: 0-1
            filter out divergence lower than cutoff percentile
        to_cpdb: bool
            whether output as CellphoneDB microenvironment file

        Returns
        ---------
        output microenvironment file

        """
    # init
    microenv = pd.Series(df_adjacency.columns, index=["Microenv_" + str(cell).replace(' ', '_')
                                                      for cell in df_adjacency.columns])
    # assign top "cutoff" percent divergence value to cutoff value
    arr_adjacency = np.array(df_adjacency).ravel()
    arr_adjacency_nonzero = []
    ## 把不是0的取出的
    for i in arr_adjacency:
        if i != 0:
            arr_adjacency_nonzero = np.append(arr_adjacency_nonzero, [i])
    ## 找到半分位点
    cutoff_v = np.percentile(np.sort(arr_adjacency_nonzero), (1-cutoff)*100)

    # find microenvironment from consensus mst
    ## 这里的microenv就是记录了不同的微环境，所谓的微环境就是
    ## 对于每一个细胞类型定义一个微环境，距离这个细胞类型进的就说明是属于这个微环境的。
    for cell in df_adjacency.columns:
        index = "Microenv_" + str(cell).replace(' ', '_')
        # find non-zero element and correspondent cell type
				## 找出每一列的值大于0且大于cutoff的，然后把对应的细胞类型记录下来
        non_zero_index = df_adjacency[cell].loc[
            (df_adjacency[cell] != 0) & (df_adjacency[cell] > cutoff_v)
            ].index.values
        # add interacting cell type
        if len(non_zero_index) != 0:
            microenv[index] = np.append(cell, non_zero_index)
        else:
            microenv.drop(index, inplace=True)

    ## 将结果转成cellphonedb的input文件
    if to_cpdb:
        out_csv_df = pd.DataFrame(columns=['cell_type', 'microenvironment'])
        for k, v in microenv.items():
            v = v.reshape(len(v), 1)
            k = np.array(str(k)).repeat(len(v)).reshape(len(v), 1)
            out_csv_df = pd.concat(
                [out_csv_df, pd.DataFrame(np.hstack([v, k]), columns=['cell_type', 'microenvironment'])],
                axis=0, ignore_index=True)
        micro_output_path = os.path.join(out_path, 'out', 'table','microenvironment.csv')
        out_csv_df.to_csv(path_or_buf=micro_output_path, header=True, index=False, sep=',')
        return out_csv_df
    else:
        return microenv

## 根据平均的js距离做热图
def divergence_clustermap(matrix, out_path):
    """ plot divergence clusterd heatmap

    Parameters
    ---------
    matrix: dataframe / matrix
        divergence matrix
    out_path: string
        output directory

    Returns
    ---------
    show and save picture

    """
    ## 计算了-log
    matrix_log = -np.log2(matrix)
    # transfre -log0 (infinite) to the max value of column
    matrix_log = matrix_log.replace(np.inf, np.nan)
    matrix_log = matrix_log.replace(np.nan, matrix_log.stack().max())
    # save csv and json
    output2echarts_divergence_cluster(matrix_log, out_path=out_path)
    # draw
    sns_plot = sns.clustermap(
        matrix_log,
        cmap="YlOrRd",
        cbar_pos=[.8, .55, .02, .2],
        dendrogram_ratio=0.1,
        method="ward")
    # mask upper triangle
    # mask = np.triu(np.ones_like(matrix_log))
    # values = sns_plot.ax_heatmap.collections[0].get_array().reshape(matrix_log.shape)
    # new_values = np.ma.array(values, mask=mask)
    # sns_plot.ax_heatmap.collections[0].set_array(new_values)
    # set left dendrogram invisible
    sns_plot.ax_row_dendrogram.set_visible(False)
    # set y axis ticks left
    sns_plot.ax_heatmap.yaxis.set_ticks_position("left")
    plt.show()
    # save figure
    clustermap_fig_output_path = os.path.join(out_path, 'out', 'pdf', 'cell_types_JSD.pdf')
    sns_plot.savefig(clustermap_fig_output_path)

def output2echarts_divergence_cluster(matrix, out_path):
    """ plot divergence clusterd heatmap

    Parameters
    ---------
    matrix: dataframe / matrix
        divergence matrix
    out_path: string
        output directory


    Returns
    ---------
    cluster matrix

    """
    # precompute linkage matrix
    Z = ward(matrix)
    index_begin = len(matrix.index)
    rearrange = {}
    for i in range(Z.shape[0]):
        x = [Z[i, 0]] if Z[i, 0] <= Z.shape[0] else rearrange[int(Z[i, 0])]
        y = [Z[i, 1]] if Z[i, 1] <= Z.shape[0] else rearrange[int(Z[i, 1])]
        x.extend(y)
        rearrange[index_begin] = x
        index_begin += 1
    cluster_list = rearrange[list(rearrange.keys())[-1]]
    cluster_matrix = matrix.iloc[cluster_list, cluster_list]
    cluster_matrix.index = cluster_matrix.columns
    # save as csv
    clustermap_csv_output_path = os.path.join(out_path, 'out', 'table','cell_types_JSD.csv')
    cluster_matrix.to_csv(path_or_buf=clustermap_csv_output_path, header=True, index=False, sep=',')
    # save as json
    cell_types = list(cluster_matrix.columns)
    log_jsd = []
    for i in range(len(cell_types)):
        for j in range(len(cell_types)):
            log_jsd.append([i, j, cluster_matrix.iloc[i, j]])
    out_dict = {
        "cell_types": cell_types,
        "log_jsd": log_jsd
    }
    clustermap_json_output_path = os.path.join(out_path, 'out', 'json','cell_types_JSD.json')
    file = open(clustermap_json_output_path, "w")
    json.dump(out_dict, file)
    file.close()

def network_draw(df_adjacency, out_path, node_size=20, edge_width=0.5):
    """ plot network graph

    Parameters
    ---------
    df_adjacency: dataframe / matrix
        adjacency matrix
    name: string
        picture name
    node_size: control node size
    edge_width: control edge width
    out_path: output pathway

    Returns
    ---------
    show and save picture

    """
    # create
    G = nx.from_pandas_adjacency(df_adjacency)
    # define position of nodes and labels
    pos = nx.random_layout(G)
    labels_pos = {}
    for key, value in pos.items():
        labels_pos[key] = (value[0], value[1])

    plt.figure(figsize=(6, 6))
    plt.axis('off')
    # save as csv
    network_csv_output_path = os.path.join(out_path, 'out', 'table','cell_types_mst_network.csv')
    df_adjacency.to_csv(network_csv_output_path)
    # save as json
    output2echarts_dict_graph(G, out_path)
    # plot nodes and edges of network graph
    nx.draw_networkx_nodes(G, pos=pos,
                           node_size=[1 * node_size * (item[1] + 1) for item in G.degree()],
                           label=True, node_color=np.array(plt.cm.YlOrRd(0.6)).reshape(1,-1))
    nx.draw_networkx_edges(G, pos=pos,
                           edge_color=[d["weight"] for (u, v, d) in G.edges(data=True)],
                           width=[d["weight"] * edge_width for (u, v, d) in G.edges(data=True)],
                           edge_cmap=plt.cm.YlOrRd)
    nx.draw_networkx_labels(G, pos=labels_pos, font_size=6, font_weight="bold")
    # plt.show(block=False)

    # save figure
    network_fig_output_path = os.path.join(out_path, 'out', 'pdf','cell_types_mst_network.pdf')
    plt.savefig(network_fig_output_path, pad_inches=0.1)
    plt.show()

def output2echarts_dict_graph(G, out_path):
    nodes = []
    edges = []
    for i in list(G.degree()):
        nodes.append({"id": i[0], "name": i[0], "degree": i[1]})
    for i in list(G.edges(data=True)):
        edges.append({"source": i[0], "target": i[1], "weight": i[2]["weight"]})
    out_dict = {
        "nodes": nodes,
        "edges": edges
    }
    network_json_output_path = os.path.join(out_path, 'out', 'json','cell_types_mst_network.json')
    file = open(network_json_output_path, "w")
    json.dump(out_dict, file)
    file.close()
