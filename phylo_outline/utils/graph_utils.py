from networkx import Graph as GraphNX
import networkx as nx
from networkx.drawing.nx_pylab import draw
from networkx_viewer import Viewer
from phylo_outline.outline.graph import Graph as MyGraph
import matplotlib.pyplot as plt


def networkx(graph: MyGraph) -> GraphNX:
    nx_graph = GraphNX()

    pos = {}
    labels = {}
    for v in graph.nodes():
        nx_graph.add_node(v.id())
        pos[v.id()] = v.location
        if v.label is not None:
            labels[v.id()] = v.label

    for e in graph.edges():
        nx_graph.add_edge(e.src().id(),e.tar().id())

    if True:
        plt.figure()
        plt.axes().axis('off')

        nx.drawing.nx_pylab.draw_networkx_edges(nx_graph, pos, width=1)
        nx.drawing.nx_pylab.draw_networkx_nodes(nx_graph, pos, node_size=1, node_color='w',alpha=0.4, node_shape='d')
        nx.drawing.nx_pylab.draw_networkx_labels(nx_graph, pos, labels=labels, font_size=12, font_family='sans-serif')
        # draw(nx_graph,pos,labels=labels)
        plt.show()
    else:
        app = Viewer(nx_graph)
        app.mainloop()


    return nx_graph



