import sys
import networkx as nx
import matplotlib.pyplot as plt

amino_acids = {'C','A','G','C','V','T','I','L','M','F','Y','W','H','K','S','P','N','D','E','R','Q'}
hydrophobic = {'C','A','G','C','V','T','I','L','M','F','Y','W','H','K'}
aromatic = {'F','Y','W','H'}
aliphatic = {'I','L','V'}
tiny = {'A','G','C','S'}
small = {'A','G','C','S','C','V','P','T','N','D'}
polar = {'K','H','R','D','E','Q','W','Y','T','C','S','N'}
charged = {'K','H','R','D','E'}
charged_pos = {'H','K','R'}
charged_neg = {'D','E'}


#precompute a dictionary of residue mappings
#this rule generator allows a single score for an exact match
#and differential scores for mismatches based on the residue membership is 
#a provided set, i.e. check if residues are both hydrophobic
def rule_generator(full_set, conditional_set, score_match,score_mismatch_both_match_condition, score_mismatch_one_matches_condition, score_mismatch_none_match_condition ):
    rule_dictionary = {}
    for r1 in full_set:
        for r2 in full_set:
            key = r1 + r2            
            #this rule generator writes each pairwise score in both directions
            #so avoid double checking
            if not (r1,r2) in rule_dictionary:
                if r1==r2:
                    rule_dictionary[key] = score_match
                else:
                    r1_in_conditional_set = r1 in conditional_set
                    r2_in_conditional_set = r2 in conditional_set
                    if r1_in_conditional_set and r2_in_conditional_set:
                        score = score_mismatch_both_match_condition
                    elif r1_in_conditional_set or r2_in_conditional_set:
                        score = score_mismatch_one_matches_condition
                    else:
                        score = score_mismatch_none_match_condition
                    rule_dictionary[key] = score
                    rule_dictionary[r2+r1] = score
    return rule_dictionary
    
matching_rule = rule_generator(amino_acids, hydrophobic,5,1,-5,0)

#function to score residue pairs between sequences
#the body of the function is just a dictionary lookup so 
#the rule may be updated by changing the dictionary
def pairwise_scoring_function(residue1,residue2):
    return matching_rule[residue1+residue2]
        
def compute_best_scored_path_for_each_node(g,m,n):
    g.node[(0,0)]["best_score"]=0    
    for j in range(0,m):
        for k in range(0,n):
            if j>0 or k>0:
                best_score = float("-inf")
                best_predecessors =[]
                for source_node_key,edge_dict in g.pred[(j,k)].items():
                    incoming_path_value = edge_dict["weight"] + g.node[source_node_key]["best_score"]
                    if incoming_path_value >  best_score:
                        best_score = incoming_path_value
                        best_predecessors = [source_node_key]
                    elif incoming_path_value == best_score:
                        #there can be more than one best path
                        best_predecessors += [source_node_key]
                current_node = g.node[j,k]
                current_node["best_score"]=best_score
                current_node["best_predecessors"]=best_predecessors
                
def create_alignment_digraph(seq1,seq2):
    m = len(seq1)+1
    n = len(seq2)+1
    g = nx.DiGraph()
    g.add_nodes_from([ (j,k) for j in range(0,m) for k in range(0,n)])
    
    indel_weight = -2
    del_edges = [((j,k),(j+1,k),indel_weight) for j in range(0,m-1) for k in range(0,n)]
    g.add_weighted_edges_from(del_edges)
    
    ins_edges = [((j,k), (j,k+1),indel_weight) for j in range(0,m) for k in range(0,n-1)]
    g.add_weighted_edges_from(ins_edges)
    
    align_edges =[((j,k),(j+1,k+1), pairwise_scoring_function(seq1[j],seq2[k])) for j in range(0,m-1) for k in range(0,n-1)]

    g.add_weighted_edges_from(align_edges)
    
    compute_best_scored_path_for_each_node(g,m,n)
                
    return g


#generate coordinates to position the nodes in a grid. 
def get_positions(m,n,grid_cell_size=10):
    top_y = grid_cell_size * (n-1)
    return {  (j,k) : (j*grid_cell_size,top_y-k*grid_cell_size) for j in range(0,m) for k in range(0,n)}


#just fetches a big list of any edges that are in a best path, does not preserve the individual paths
#this result is just used for highlighting edges graphically, not for aligning the sequences
def get_edges_in_best_paths(g,into_node_key):
    sink = g.node[into_node_key]
    if 'best_predecessors' in sink:
        predecessors = sink['best_predecessors']
        edges = []
        for pred in predecessors:
            edges = edges + [(pred,into_node_key)] + get_edges_in_best_paths(g,pred)
        return edges    
    else:
        return []

#this fetches a list of any and all best paths
def get_best_paths(g,into_node_key):
    sink = g.node[into_node_key]
    if 'best_predecessors' in sink:
        predecessors = sink['best_predecessors']
        paths = []
        for pred in predecessors:
            parental_paths = get_best_paths(g,pred)
            if len(parental_paths)>0:
                for p in parental_paths:
                    new_entry = [into_node_key] + p
                    paths.append(new_entry) #new_entry)

            else:
                paths.append(pred)
        return paths    
    else:
        return [[into_node_key]]

#draw a small arrows in each cell indicating the direction from which
#the best path(s) came
def draw_traceback_indices(g,positions,seq1,seq2, cell_size):
    hfont = {'fontname':'Wingdings 3','fontsize':'14'}
    for j in range(0,len(seq1)+1):
        for k in range(0,len(seq2)+1):
            pos = positions[(j,k)]
            node = g.node[(j,k)]
            if 'best_predecessors' in node:
                pred = node['best_predecessors']
                if len(pred)>0:
                    arrows = ""
                    for p in pred:
                        if p[0]==j and p[1]==k-1:
                            arrows+="h" #in Wingdings 3 h is up arrow
                        if p[0]==j-1 and p[1]==k:
                            arrows+="f" #in Wingdings 3 f is left arrow
                        if p[0]==j-1 and p[1]==k-1:
                            arrows+="j" #in Wingdings 3 f is NW arrow
                    if len(arrows)>0:
                        plt.text(pos[0] - cell_size * 0.45 , pos[1] +  cell_size * 0.65,s=arrows,**hfont,color="m") 
        

def draw_alignment_grid(g,seq1,seq2,filename="alignment.png",plot=False):
    cell_size = 10
    positions = get_positions(len(seq1)+1,len(seq2)+1,cell_size)
    #edge_positions = {key:(positions[key][0]-20, positions[key][1]+20) for key in positions}
    edge_positions = positions

    plt.rcParams['figure.figsize'] = (12.0, 10.0)

    plt.axis('off')

    labels= {key:data['best_score'] for key,data in g.nodes(data=True) }
    nx.draw_networkx(g,pos=positions,labels=labels,node_color="lightblue")

    labels = nx.get_edge_attributes(g,'weight')

    nx.draw_networkx_edge_labels(g, pos=edge_positions,edge_labels=labels,font_size=7)

    draw_traceback_indices(g,positions,seq1,seq2, cell_size)

    origin_position = positions[(0,0)]
    for r in range(0,len(seq1)):
        plt.text(origin_position[0] + r * cell_size + cell_size/2,origin_position[1] +cell_size,s=seq1[r], bbox=dict(facecolor='red', alpha=0.5),horizontalalignment='center')
        
    for r in range(0,len(seq2)):
        plt.text(origin_position[0] - cell_size,origin_position[1]- r*cell_size - cell_size/2,s=seq2[r], bbox=dict(facecolor='red', alpha=0.5),horizontalalignment='center')
        
    #sink = g.node[m-1,n-1]        

    best_path_edges = get_edges_in_best_paths(g,(len(seq1),len(seq2)))


    nx.draw_networkx_edges(g,pos=positions,edgelist=best_path_edges,edge_color="r",width=8,alpha=0.5)
    plt.savefig(filename)
    if plot:
        plt.show() 


def map_alignments(g,seq1,seq2):
    best_paths = get_best_paths(g,(len(seq1),len(seq2)))
    
    alignments = []
    for path in best_paths:

        aligned_seq1 = ""
        aligned_seq2 = ""
        path = list(reversed(path))
        for i in range(0,len(path)-1):
            source = path[i]
            sink = path[i+1]
            move_right = sink[0]-source[0]==1
            move_down = sink[1]-source[1]==1
            if move_right:
                if move_down:  #the diagonal case, both residues are aligned
                    aligned_seq1 += seq1[source[0]]
                    aligned_seq2 += seq2[source[1]]
                else:
                    aligned_seq1 += seq1[source[0]]
                    aligned_seq2 += "_"
            elif move_down:
                aligned_seq1 += "_"
                aligned_seq2 += seq2[source[1]]
        alignments += [(aligned_seq1,aligned_seq2)]
    return alignments

def display_alignments(g,seq1,seq2):
    r = map_alignments(g,seq1,seq2)
    
    l = len(r)
    if l>1:
        plural="s"
    else:
        plural=""
        
    print("%d best path%s found:" % (l,plural))
    print()
    for a in r:
        print(a[0])
        print(a[1])
        print()


def align_and_display(seq1,seq2,filename="",plot=False):
    g = create_alignment_digraph(seq1, seq2)    

    draw_alignment_grid(g,seq1,seq2, filename,plot)            

    display_alignments(g,seq1,seq2)