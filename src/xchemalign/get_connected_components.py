import networkx as nx
import numpy as np
from loguru import logger


def get_connected_components(arr):

    g = nx.Graph()
    for idx, conn in np.ndenumerate(arr):
        x, y = idx
        if conn:
            g.add_edge(x, y)

    cliques = list(nx.enumerate_all_cliques(g))
    logger.debug(cliques)

    logger.debug(g)
