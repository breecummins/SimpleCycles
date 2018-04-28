import DSGRN
import networkx as NX
import sys, time, itertools, operator

def separateMultipleOrders(names,extrema):
    ''' 
    If there are two or more maxes in a row, then delete every one but one.
    Do the same for mins and multiple groups of greater than 1 max or min.
    This returns a list of sequences with alternating maxes and mins in each variable. 

    '''
    def group_inds(l):
        it  = itertools.groupby(l,operator.itemgetter(1))
        return [ (key,[item[0] for item in subiter ]) for key, subiter in it ]

    extrema_inds = []
    for name in names:
        # get all the extrema for variable name
        nf = filter(lambda x: x[1].startswith(name),enumerate(extrema))
        # record the index of each extremum and whether or not it's a max or min
        ext = [(e[0],e[1].split()[1]) for e in nf]
        # group by consecutive maxes or consecutive mins and record the bad sequences
        grouped_inds = group_inds(ext)
        ginds = [ i[1] for i in grouped_inds]
        # make sure that consecutive blocks wrap around the sequence
        if grouped_inds and (grouped_inds[0][0] == grouped_inds[-1][0]):
            ginds = [ginds[0]+ginds[-1]] + ginds[1:-1]
        # filter for blocks with more than 1 identical extrema
        inds = filter(lambda x : len(x) > 1, ginds)
        if inds: extrema_inds.extend(inds)
    # if the sequence has alternating maxes and mins to start with, record it as good
    if not extrema_inds:
        new_extrema_list = [extrema]
    # if the sequence has consecutive blocks of identical extrema, find every compatible 
    # sequence with alternating maxes and mins in each variable
    else:
        new_extrema_list = []
        X = [ x for ext in extrema_inds for x in ext ]
        for prod in itertools.product(*extrema_inds):
            #remove extra indices NOT IN prod
            new = tuple([v for i,v in enumerate(extrema) if (i not in X) or (i in prod)])
            if new not in new_extrema_list: new_extrema_list.append(new)                 
    return new_extrema_list

def orderedExtrema(names,labeled_cycles):
    # need to write function based on paths of the form -m-, ---, M--, etc.
    allextrema = set([])
    for cyc in labeled_cycles:
        extrema = []
        for c in cyc:
            if "*" in c: raise ValueError("Debug: * in label.")
            elif "M" in c: extrema.append(names[c.index("M")]+" max")
            elif "m" in c: extrema.append(names[c.index("m")]+" min")
            else: pass  
        new_extrema_list = separateMultipleOrders(names,tuple(extrema))
        for e in new_extrema_list: allextrema.add( e )
    return allextrema

def findCycles(digraph):
    # graph is nx.DiGraph object
    cycles = NX.simple_cycles(digraph)
    cycles = [cyc+[cyc[0]] for cyc in cycles] #first element is left off of the end in simplecycles() output
    labeled_cycles = set([tuple([digraph.edge[u][v]["label"] for (u,v) in zip(cyc[:-1],cyc[1:])]) for cyc in cycles])
    return labeled_cycles

def makeNXDigraph(domaingraph,nodes=None,edges=None):
    ''' 
    Make networkx digraph in order to use the networkx library to find simple cycles.

    '''
    if nodes is None:
        # get nodes
        nodes = range(domaingraph.digraph().size())
    if edges is None:
        # get edges
        edges=[ (i,a) for i in nodes for a in domaingraph.digraph().adjacencies(i) ]
    # attach labels to edges
    searchgraph = DSGRN.SearchGraph(domaingraph)
    MR = DSGRN.MatchingRelation(domaingraph.dimension())
    edgelabels = { (i,a) : MR.edge_labelstring(searchgraph.event(i,a)) for (i,a) in edges }
    # add nodes and edges to digraph
    G = NX.DiGraph()
    G.add_nodes_from(nodes)
    for edge,label in edgelabels.items(): G.add_edge(edge[0],edge[1],label=label)
    return G


# def findAllOrderedExtrema(networkfile=None,networkspec=None):
#     if networkfile:
#         network = DSGRN.Network(networkfile)
#     elif networkspec:
#         network = DSGRN.Network()
#         network.assign(networkspec)
#     else:
#         raise ValueError("No input network.")
#     names = [network.name(i) for i in range(network.size())]
#     paramgraph = DSGRN.ParameterGraph(network)
#     paths = set([])
#     mastercycles = set([])
#     for paramind in range(paramgraph.size()):
#         domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
#         cycles = findCycles(makeNXDigraph(domaingraph))
#         for cyc in cycles:
#             # check for identical cycles at different starting nodes
#             if cyc not in mastercycles:
#                 mastercycles.update(makeMasterCycles(cyc))
#                 paths.update(orderedExtrema(names,cycles))  
#     return set(paths)

def notInCyclicPermutations(x,cycle):
    return x not in [ tuple(list(cycle[n:]) + list(cycle[:n])) for n in range(len(cycle)) ]

def removeCyclicPermutations(extrema,paths):
    new_paths = paths.copy()
    for e in extrema:
        # check if any existing paths have the same length and set value
        same_len = [ p for p in new_paths if len(e)==len(p) and set(e)==set(p) ]
        # if not, then add to path list
        if not same_len:
            new_paths.add( e )
        # if so, then check for cyclic permutations
        else:
            different = True
            while different and same_len: different = notInCyclicPermutations( e, same_len.pop() )
            if different: new_paths.add( e )
    return new_paths


def findAllOrderedExtrema_Morsesets(networkfile=None,networkspec=None):
    if networkfile:
        network = DSGRN.Network(networkfile)
    elif networkspec:
        network = DSGRN.Network()
        network.assign(networkspec)
    else:
        raise ValueError("No input network.")
    names = [network.name(i) for i in range(network.size())]
    paramgraph = DSGRN.ParameterGraph(network)
    paths = set([])
    start = time.time()
    for paramind in range(paramgraph.size()):
        if time.time()-start >= 2:
            print("{} / {} parameters analyzed\n".format(paramind,paramgraph.size()))
            start = time.time()
        domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
        morsedecomposition = DSGRN.MorseDecomposition(domaingraph.digraph())
        morsegraph = DSGRN.MorseGraph()
        morsegraph.assign(domaingraph,morsedecomposition)
        poset = morsedecomposition.poset()
        for i in range(0,morsedecomposition.poset().size()):
            ms = morsedecomposition.morseset(i)
            if len(ms) > 1 and morsegraph.annotation(i)[0] == "FC" and len(poset.children(i)) == 0:
                morseedges = [ (j,a) for j in ms for a in domaingraph.digraph().adjacencies(j) if a in ms ]
                digraph = makeNXDigraph(domaingraph,ms,morseedges)
                print("Nodes: {}".format(digraph.number_of_nodes()))
                print("Edges: {}".format(digraph.size()))
                cycles = findCycles(digraph)
                print("Have cycles.")
                k = 0
                for c in cycles: k+=1
                print("Number cycles: {}".format(k))
                sys.exit()
                # debugging try-except block
                try: 
                    C = max(len(c)-1 for c in cycles)
                    if C > len(ms):
                        print("morse set: {}, max cycle: {}".format(len(ms),C))
                        raise ValueError("Nodes in cycle exceeds nodes in Morse set.")
                except: pass
                extrema  = orderedExtrema(names,cycles)
                paths = removeCyclicPermutations(extrema,paths)
                print(paths)
                sys.stdout.flush()
    return set(paths)

def findAllOrderedExtremaDomainGraph(paramlist=None,networkfile=None,networkspec=None):
    if networkfile:
        network = DSGRN.Network(networkfile)
    elif networkspec:
        network = DSGRN.Network()
        network.assign(networkspec)
    else:
        raise ValueError("No input network.")
    names = [network.name(i) for i in range(network.size())]
    paramgraph = DSGRN.ParameterGraph(network)
    paths = []
    start = time.time()
    if not paramlist:
        paramlist=range(paramgraph.size())
    for paramind in paramlist:
        if time.time()-start >= 2:
            print("{} / {} parameters analyzed\n".format(paramind,paramgraph.size()))
            start = time.time()
        domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
        digraph = makeNXDigraph(domaingraph)
        cycles = findCycles(digraph)
        extrema  = orderedExtrema(names,cycles)
        p = removeCyclicPermutations(extrema,set([]))
        paths.append(list(p))
    return paths

def test_multiple_extrema():
    names = ['x','y','z']
    extrema = ('x max','y max','z min','y max','x min','y min','z max')
    print(set(separateMultipleOrders(names,extrema)) == set([('x max','z min','y max','x min','y min','z max'),('x max','y max','z min','x min','y min','z max')]))

    extrema = ('x max','y max','z min','y max','x min','x min', 'y min','z max')
    # print separateMultipleOrders(names,extrema_list)
    print(set(separateMultipleOrders(names,extrema)) == set([('x max','z min','y max','x min','y min','z max'),('x max','y max','z min','x min','y min','z max')]))


    extrema = ('x max','y max','z min','y max','x min','x min', 'y min', 'x min', 'z max')
    # print separateMultipleOrders(names,extrema_list)
    print(set(separateMultipleOrders(names,extrema)) == set([('x max','z min','y max','x min','y min','z max'),('x max','z min','y max','y min','x min','z max'),('x max','y max','z min','x min','y min','z max'),('x max','y max','z min','y min','x min','z max')]))


if __name__ == "__main__":
    # netspec0 = "X : ~Z\nY : ~X\nZ : ~Y"
    # netspec1 = "x1 : ~z1 : E\ny1 : ~x1 : E\nz1 : ~y1 : E\nx2 : ~z2 : E\ny2 : (~x1)(~x2) : E\nz2 : ~y2 : E"
    # netspec2 = "x1 : ~z1 : E\ny1 : (~x1)(~x2) : E\nz1 : ~y1 : E\nx2 : ~z2 : E\ny2 : (~x1)(~x2) : E\nz2 : ~y2 : E"
    # netspec3 = "x1 : ~z1 : E\ny : (~x1)(~x2) : E\nz1 : ~y : E\nx2 : ~z2 : E\nz2 : ~y : E"
    # netspec4 = "X : (~Z)(Y) : E\nY : ~X : E\nZ : ~Y : E"

    # num = sys.argv[1]
    # ns = eval("netspec" + num)
    # print "\n" + ns + "\n\n"
    # sys.stdout.flush()

    # # for c in list(findAllOrderedExtrema_Morsesets(networkspec=ns))[:6]:
    # #     print c

    # orders = findAllOrderedExtrema_Morsesets(networkspec=ns)
    # print len(orders)
    # for o in sorted(list(orders)):
    #     print o

    # # test_multiple_extrema()

    paths = findAllOrderedExtremaDomainGraph(paramlist = [5,9,19,23,30,38,44,52],networkfile="syn_net_5D_20180425.txt")

    import json
    D = dict(zip([5,9,19,23,30,38,44,52],paths))
    json.dump(D,open("syn_net_5D_20180425_cycles.json","w"))
    

