// Last edited 24/08/20

// This is a C++ implementation of Hierarchical Density Based Spatial Clustering of Applications with Noise (HDBSCAN). The original
// version was developed by B. Harvey; a second, more lightweight version made for in-script implementation was drafted by T. Hankins.

//_________________________________________________________________________________________________
// HDBSCAN Usage
// - The clusterer is controlled by a set of parameters that can be adjusted after instantiation.
//      - k: number of nearest neighbors in the standard metric search
//          - Adjusting this can help with convergence on secondary peaks in a spectrum.
//      - minimum_cluster_size: number of data points required for a cluster to form
//      - alpha: weighting factor for choosing optimal clusters; larger alpha = preferentially more clusters
//      - allow_single_cluster: whether or not a single cluster should be allowed
//  
// - Each property can be controlled by calling set_[](), where [] is the name of the variable (e.g.,
//   set_minimum_cluster_size()).
//_________________________________________________________________________________________________

#ifndef HDBSCAN__H
#define HDBSCAN__H

#include <iostream>
#include <vector>

#include "TMath.h"
#include "TStopwatch.h"

class Point;

class Node {
    private:
        Node *parent {};
        std::vector<Node *> children {};

        double birth_lambda {}; // 1/(distance which caused cluster to separate into its own). Unless otherwise told, 
                                // assume infinitely far away from other Nodes
        double death_lambda {std::numeric_limits<double>::infinity()}; // 1/(distance which caused cluster to separate into smaller clusters)
        double stability {-1};

        int size {};
        std::vector<Point *> member_points {};

    public:
        Node() {};
        Node(Node *, Node *, double);
        Node(Point *);
        ~Node();

        // Getters
        Node *get_child(int index);
        Node *get_parent() { return parent; }
        int get_size() { return size; }
        double get_stability();
        double get_birth_lambda() { return birth_lambda; }
        double get_death_lambda() { return death_lambda; }
        int get_number_leaf_nodes();
        int get_number_leaf_data();
        std::vector<Point *> get_contained_points();
        std::vector<Point *> get_member_points() { return member_points; }
        std::vector<Point *> get_exemplar_points();

        // Setters
        void set_size(int size) {this->size = size; }
        void set_default_stability();
        void set_stability(double stability) { this->stability = stability; }
        void set_parent(Node *parent) { this->parent = parent; }
        void set_birth_lamdba(double lambda) { birth_lambda = lambda; }
        void set_death_lambda(double lambda) { death_lambda = lambda; }
        void add_child(Node *cluster) { children.push_back(cluster); }
        void set_children(Node *cluster_a, Node *cluster_b) { children = {cluster_a, cluster_b}; }
        void add_member_point(Point *point) { member_points.push_back(point); ++size; }
        
        // Other
        void find_points(std::vector<Point *> &);
        void find_exemplar_points(std::vector<Point *> &);
        Node *find_condensed();
        bool is_leaf();
};

//__________________________________________________________________________________________________________________________________

class Point {
    private:
        std::vector<double> coordinate {}; // point location
        double core_distance {}; // measure of inverse density, distance to k-th nearest neighbor

        std::vector<Point *> neighbors {}; // for reconnecting pruned forests

        bool connected {false};

        int id {-1}; // set to -1 until a clusterer says otherwise
        Node *leaf_cluster {}; // points to cluster that originally only contained this point
        double lambda {-1}; // 1/distande threshold required to add point to cluster

        std::vector<double> cluster_confidence_levels {};
        std::vector<Node *> final_containing_clusters {};

    public:
        Point() {};
        Point(std::vector<double>);
        // Point(const Point &);
        ~Point() { std::cout << "Point::~Point()" << std::endl; }

        // Getters
        double get_coordinate(int dim);
        double get_x() { return get_coordinate(0); }
        double get_y() { return get_coordinate(1); }
        double get_z() { return get_coordinate(2); }
        int get_number_dimensions() { return static_cast<int>(coordinate.size()); }
        double get_core_distance() { return core_distance; }
        double get_relative_density() { return TMath::Power(core_distance, -1.0 * get_number_dimensions()); }
        std::vector<double> get_cluster_confidence_levels() { return cluster_confidence_levels; }
        std::vector<Node *> get_final_containing_cluster() { return final_containing_clusters; }
        int get_cluster_id() { return id; }
        Node *get_leaf_cluster();
        double get_lambda() { return lambda; }
        std::vector<Point *> get_neighbors() { return neighbors; }

        double get_euclidean_distance(Point *);
        double get_mutual_reachability_distance(Point *);
        bool is_connected() { return connected; }

        std::vector<Point *> *get_neighbors_address() { return &neighbors; }

        // Setters
        void set_point(std::vector<double> coordinate);
        void set_core_distance(double core_distance) { this->core_distance = core_distance; }
        void set_cluster_id(int id) { this->id = id; }
        void set_leaf_cluster(Node *leaf_cluster) { this->leaf_cluster = leaf_cluster; }
        void set_lambda(double lambda) { this->lambda = lambda; }

        // Other
        void add_neighbor(Point *p) { neighbors.push_back(p); }
        void mark_connected() { connected = true; }
        Node *get_cluster();
};

//__________________________________________________________________________________________________________________________________

class Edge {
    private:
        Point *connected_spatial_points[2] = {0, 0};
        double length {};

    public:
        Edge();
        ~Edge() {std::cout << "Edge::~Edge()" << std::endl; }
        Edge(Point *, Point *);
        Edge(Point *, Point *, double);

        double get_length() { return length; }
        bool contains_point(Point *p) { return (connected_spatial_points[0] == p || connected_spatial_points[1] == p); }
        Point *get_point(int index) {
            if(index == 0 || index == 1) {
                return connected_spatial_points[index];
            } else {
                return connected_spatial_points[0];
            }
        }
};

//__________________________________________________________________________________________________________________________________

class Vertex {
    private:
        Point *datum {};
        Vertex *representative {};
        Vertex *merging {};
        double delta {std::numeric_limits<double>::max()};
        Vertex *foreign_neighbor {};

    public:
        Vertex() {};
        ~Vertex() {};
        Vertex(Point *p) : datum(p), representative(this) {};

        // Getters
        Point *get_datum() { return datum; }
        Vertex *get_parent() { return representative; }
        double get_delta() { return delta; }
        Vertex *get_foreign_neighbor() {return foreign_neighbor; }
        Vertex *get_merging_vertex() { return merging; }
        Vertex *get_component() {
            if(this->get_parent() != this) {
                this->set_parent(this->get_parent()->get_component());
                return this->get_parent();
            } else {
                return this;
            }
        }

        // Setters
        void set_datum(Point *p) { datum = p; }
        void set_parent(Vertex *v) { representative = v; }
        void set_delta(double value) { delta = value; }
        void set_foreign_neighbor(Vertex *v) { foreign_neighbor = v; }
        void set_merging_vertex(Vertex *v) { merging = v; }

};

//__________________________________________________________________________________________________________________________________

class Ball {
    private:
        double radius {};
        int number_of_points {};
        std::vector<Vertex *> contained_points {};
        Point *pivot {};
        Ball *child_1{};
        Ball *child_2 {};
        bool is_mrd_tree {};
        TString distance_function {"Euclidean"};
        double distance_to_component {};
        Vertex *contained_component {};

        // Helpers
        double get_standard_distance(Point *, Point *);
        double get_mutual_reachability_distance(Point *a, Point *b) {
            return std::max(std::max(a->get_core_distance(), b->get_core_distance()), get_standard_distance(a, b));
        }
        void initialize_standard_ball(std::vector<Vertex *>, TString, int);
        void initialize_mrd_ball(std::vector<Vertex *>, TString, int, Ball *, int);
        void find_contained_points(std::vector<Vertex *> &);

    public:
        Ball() {};
        ~Ball() {};

        // Standard metric
        Ball(std::vector<Point *>, TString, int);
        Ball(std::vector<Vertex *>, TString, int);

        // Mutual reachability depending on a standard metric
        Ball(std::vector<Point *>, TString, int, Ball *, int);
        Ball(std::vector<Vertex *>, TString, int, Ball *, int);

        // Getters
        double get_distance_to_component() { return distance_to_component; }
        Vertex *get_contained_component() {
            if(contained_component) return contained_component->get_component();
            return 0;
        }
        Ball *get_child_1() { return child_1; }
        Ball *get_child_2() { return child_2; }
        int get_number_of_points() {return number_of_points; }
        Point *get_pivot() { return pivot; }
        double get_radius() { return radius; }
        std::vector<Vertex *> get_contained_data();
        double get_distance(Point *, Point *);
        double get_distance(Vertex *, Vertex *);

        // Setters
        void set_metric(TString);
        void set_distance_to_component(double distance) { distance_to_component = distance; }
        void set_component(Vertex *v) { contained_component = v; }
        void set_contained_components();

        // Finders
        void find_knn(Point *, int, std::vector<Point *>&);
        void find_knn(Point *, int, std::vector<Vertex *>&);
        void find_radial_nn(Point *, double, std::vector<Point *>&);
        void find_radial_nn(Point *, double, std::vector<Vertex *>&);
        void find_radial_unsorted_nn(Point *, double, std::vector<Point *>&);
        void find_radial_unsorted_nn(Point *, double, std::vector<Vertex *>&);
        void find_foreign_nn(Vertex *, Vertex *);

        // Other
        bool is_leaf() { return this == child_1; }
        void reset_distance_to_components();
        
};

//__________________________________________________________________________________________________________________________________

class MinimumSpanningTree {
    private:
        std::vector<Edge *> edge_set {};
        std::vector<Vertex *> vertex_set {};
        bool edge_sorted {true}; // Flag to restrict sorting to only after operations which
                                 // would result in the edge set losing its sort

        // Nearest foreign neighbor querty
        void find_nearest_foreign_vertex(Vertex *, Ball *, Vertex *);
        // For dual tree Boruvka
        bool completely_connected(Ball *, Ball *);
        void add_edges(std::vector<Vertex *> &);

    public:
        MinimumSpanningTree() {}
        MinimumSpanningTree(Ball *);
        MinimumSpanningTree(std::vector<Vertex *>, std::vector<Edge *>);
        ~MinimumSpanningTree() {};

        // Getters
        int get_size() { return static_cast<int>(edge_set.size()); }
        int get_order() { return static_cast<int>(vertex_set.size()); }
        std::vector<Edge *> get_edge_set() { return edge_set; }
        Edge *get_edge_at(int i) { return edge_set.at(i); }
        std::vector<Vertex *> get_vertex_set() { return vertex_set; }
        Vertex *get_vertex_at(int i) { return vertex_set.at(i); }
        Edge *get_nth_longest_edge(int);
        bool contains_edge(Edge *);

        // Other
        void add_edge(Edge *e) {
            edge_set.push_back(e);
            edge_sorted = false;
        }
        void add_vertex(Vertex *v) { vertex_set.push_back(v); }
        void sort_edges();

        // for dual tree Boruvka
        void find_component_neighbors(Ball *, Ball *);
};

//__________________________________________________________________________________________________________________________________

class HDBSCAN {
    private:
        std::vector<Point *> data {};

        // Clusterer Parameters
        int k {};                                   // Number of nearest neighbors in the standard metric search
        int minimum_cluster_size {};                // Number of data points required for a cluster to form in the condensed tree
        double alpha {};                            // Weighting factor for choosing optimal clusters; larger alpha = preferentially more clusters;
                                                        // does not affect the calculated cluster node tree, only which clusters are selected from it
        int leaf_size {};                           // The maximum number of spatial points stored in a leaf node of a ball. This only affects
                                                        // performance. Smaller leaves will generally result in more time efficient and storage
                                                        // inefficient structures
        bool allow_single_cluster {};               // Whether or not the clustering should allow a single cluster
        TString cluster_selection_method {"EOM"};   // Determines which option to cluster with: "EOM" (default) or "Leaf"

        // Scaling factors
        std::vector<double> scales {};              // Determines the scale of the data. Clustering data which is set on similar scales is generally
                                                        // better. Large scales weight the importance of that axis higher than other axes. Default is
                                                        // ones for each dimension. Zeros will prevent the data from being scaled.
        
        // Metric
        TString metric {"E"};                       // The metric used to determine distance. "Euclidean"/"E" is the default (and only) option.

        // Ball trees
        Ball *standard_metric_ball_root {};         // Spatial indexing tree for nearest neighbor queries
        Ball *mutual_reachability_ball_root {};     // Spatial indexing tree for nearest neighbor queres using mutual reachability distance (MRD)

        // Minimum spanning tree
        MinimumSpanningTree *minimum_spanning_tree {}; // The minimum spanning tree using MRD to connect all points

        // Cluster dendrogram
        Node *root_cluster {};                  // Head of the cluster node tree
        std::vector<Node *> leaf_clusters {};   // Leaf node clusters (have no children)
        std::vector<Node *> final_clusters {};  // The set clusters at the end of the fitting routine

        // Keep track of the minimum and maximum values for the original data so that space transformations can be performed on new data
        std::vector<double> minimum_values {};
        std::vector<double> maximum_values {};

        std::vector<TString> variable_names {};

        bool scaled {};                         // Keeps track of whether the data has been scaled yet

        void set_extrema_values();

        // Data members are deleted upon relevant parameter changes
        void dispose_mst();
        void dispose_ball_trees();
        void dispose_cluster_node_tree();

        // Routines for fitting
        void initialize_ball_trees();           // Nearest neighbor query
        void construct_minimum_spanning_tree(); // Minimum spanning tree search (uses k)
        void construct_condensed_cluster_hierarchy(); // Create dendrogram internally
        // Add edge from MST to connect clusters
        void merge_condensed_clusters(Point *, Point *, double);
        // Excess of mass cluster selection
        std::vector<Node *> optimize_clusters_by_eom(Node *);
        std::vector<Node *> get_leaf_clusters(Node *);
        void assign_clusters();                 // Add clusters to the clusters attribute by optimizing clusters
        void label_clusters();                  // Use index in clusters to assign labels to each point in data

    public:
        HDBSCAN() {}
        HDBSCAN(std::vector<Point *> &data) : HDBSCAN(data, "") {}
        HDBSCAN(std::vector<Point *> &data, TString parameter_file_name);

        ~HDBSCAN() {}

        // Parameter setting/getting
        void set_parameters(); // Set default parameters
        void set_parameters(TString);

        void set_metric(TString);
        void set_k(int);
        void set_ball_tree_leaf_size(double);
        void set_minimum_cluster_size(int);
        void set_alpha(double);
        void set_cluster_selection_method(TString);
        void set_allow_single_cluster(bool allow) { allow_single_cluster = allow; }
        void set_variable_names(std::vector<TString> variable_names) { this->variable_names = variable_names; }
        void set_scales(std::vector<double>);
        void scale_data();
        void unscale_data();
        void print_parameters();

        int get_k() { return k; }
        bool get_scaled() { return scaled; }
        std::vector<double> get_scaled_coordinates(Point *);
        int get_minimum_cluster_size() { return minimum_cluster_size; }
        double get_alpha() { return alpha; }
        Ball *get_standard_metric_ball_tree() { return standard_metric_ball_root; }
        Ball *get_mutual_reachability_ball_tree() { return mutual_reachability_ball_root; }
        int get_number_points() { return static_cast<int>(data.size()); }
        std::vector<Point *> *get_data_location() { return &data; }
        TString get_variable_name(int index) { 
            return (index < static_cast<int>(variable_names.size())) ? variable_names.at(index) : Form("Variable %i", index); 
        }

        // Call the fitting routines in the proper order
        void fit();

        // Extra functionality
        std::vector<double> get_unscaled_coordinates(Point *);
        std::vector<Node *> get_final_clusters() { return final_clusters; }
        std::vector<Point *> *get_data() { return &data; }
};

//__________________________________________________________________________________________________________________________________

Node::Node(Point *p) {
    parent = this;
    children = {this, this};
    member_points = {p};
    size = 1;
}

Node::Node(Node *cluster_a, Node *cluster_b, double length) {
    parent = this;
    children = {cluster_a, cluster_b};
    death_lambda = 1.0 / length; // length is length of edge connecting Nodes x and y
    if(cluster_a->get_size() < cluster_b->get_size()) std::swap(cluster_a, cluster_b);
    size = cluster_a->get_size() + cluster_b->get_size();
}

Node::~Node() {
    if(parent != this) {
        if(parent->get_child(0) == this) {
            parent->set_children(parent, parent->get_child(1));
        }

        if(parent->get_child(1) == this) {
            parent->set_children(parent->get_child(0), parent);
        }
    }

    if(get_child(0) != this) get_child(0)->set_parent(get_child(0));
    if(get_child(1) != this) get_child(1)->set_parent(get_child(1));

    for(auto *p : member_points) p->set_leaf_cluster(0);

    parent = 0;
    children = {};
    member_points = {};
}

Node *Node::get_child(int index) {
    if(index >= static_cast<int>(children.size())) {
        std::cout << "Node::get_child() : index out of bounds, returning last element" << std::endl;
        return children.back();
    } else {
        return children.at(index);
    }
}

double Node::get_stability() {
    if(stability < 0) set_default_stability();
    return stability;
}

int Node::get_number_leaf_nodes() {
    if(is_leaf()) return 1;
    return (children.at(0)->get_number_leaf_nodes() + children.at(1)->get_number_leaf_nodes());
}

int Node::get_number_leaf_data() {
    if(is_leaf()) return size;
    return (children.at(0)->get_number_leaf_data() + children.at(1)->get_number_leaf_data());
}

std::vector<Point *> Node::get_contained_points() {
    std::vector<Point *> members {};
    find_points(members);
    return members;
}

std::vector<Point *> Node::get_exemplar_points() {
    std::vector<Point *> members {};
    find_exemplar_points(members);
    return members;
}

void Node::set_default_stability() {
    stability = 0;
    std::vector<Point *> points {};
    find_points(points);
    for(auto *p : points) stability += (std::min(p->get_lambda(), death_lambda) - birth_lambda);
}

void Node::find_points(std::vector<Point *> &points) {
    for(auto *p : member_points) points.push_back(p);

    if(get_child(0) != this) {
        get_child(0)->find_points(points);
        get_child(1)->find_points(points);
    }
}

void Node::find_exemplar_points(std::vector<Point *> &points) {
    if(this->is_leaf()) {
        for(auto *p : member_points) {
            if(p->get_lambda() >= this->get_death_lambda()) points.push_back(p);
        }
    } else {
        get_child(0)->find_exemplar_points(points);
        get_child(1)->find_exemplar_points(points);
    }
}

Node *Node::find_condensed() {
    if(get_parent() != this) {
        return get_parent()->find_condensed();
    } else {
        return this;
    }
}

bool Node::is_leaf() {
    if(static_cast<int>(children.size()) > 0) {
        return (children.at(0) == this);
    } else {
        return true;
    }
}

//__________________________________________________________________________________________________________________________________

Point::Point(std::vector<double> coordinate) {
    for(auto value : coordinate) this->coordinate.push_back(value);
}

// Point::Point(const Point &)

double Point::get_coordinate(int dim) {
    if(dim < static_cast<int>(coordinate.size())) return coordinate.at(dim);
    return 0.0;
}

Node *Point::get_leaf_cluster() {
    if(leaf_cluster) return leaf_cluster;
    return 0;
}

double Point::get_euclidean_distance(Point *p) {
    double distance_sq {};
    int number_dimensions {std::max(this->get_number_dimensions(), p->get_number_dimensions())};
    for(int i = 0; i < number_dimensions; ++i) {
        distance_sq += TMath::Power(this->get_coordinate(i) - p->get_coordinate(i), 2.0);
    }
    return TMath::Sqrt(distance_sq);
}

double Point::get_mutual_reachability_distance(Point *p) {
    return std::max(std::max(this->get_core_distance(), p->get_core_distance()), this->get_euclidean_distance(p));
}

void Point::set_point(std::vector<double> coordinate) {
    this->coordinate.clear();
    for(auto value : coordinate) this->coordinate.push_back(value);
}

Node *Point::get_cluster() {
    if(leaf_cluster) return leaf_cluster->find_condensed();
    return 0;
}

//__________________________________________________________________________________________________________________________________

Edge::Edge(Point *p1, Point *p2) {
    connected_spatial_points[0] = p1;
    connected_spatial_points[1] = p2;
    length = p1->get_mutual_reachability_distance(p2);
}

Edge::Edge(Point *p1, Point *p2, double length) {
    connected_spatial_points[0] = p1;
    connected_spatial_points[1] = p2;
    this->length = length;
}

//__________________________________________________________________________________________________________________________________

double Ball::get_standard_distance(Point *a, Point *b) {
    return a->get_euclidean_distance(b);
}

void Ball::initialize_standard_ball(std::vector<Vertex *> data, TString metric, int leaf_size) {
    number_of_points = static_cast<int>(data.size());
    distance_to_component = std::numeric_limits<double>::infinity();

    // Pick a point x in the ball, y and z will be far apart
    Vertex *x {data.back()}, *y {}, *z {};
    double max_distance {}, distance {};

    // y is the furthest point in the ball from x
    for(auto *v : data) {
        if(!y) {
            y = v;
            max_distance = get_distance(v, x);
        } else {
            distance = get_distance(v, x);
            if(distance > max_distance) {
                y = v;
                max_distance = distance;
            }
        }
    }

    max_distance = 0;
    for(auto *v : data) {
        if(v != y) {
            if(!z) {
                z = v;
                max_distance = get_distance(v, y);
            } else {
                distance = get_distance(v, y);
                if(distance > max_distance) {
                    z = v;
                    max_distance = distance;
                }
            }
        }
    }

    // Define a midpoint of y and z to be the midpoint of the ball
    std::vector<double> midpoint_coordinates {};
    if(!z) z = y; // Should only happen with one point in the ball

    for(int i = 0; i < y->get_datum()->get_number_dimensions(); ++i) {
        midpoint_coordinates.push_back((y->get_datum()->get_coordinate(i) + z->get_datum()->get_coordinate(i)) / 2.0);
    }

    pivot = new Point(midpoint_coordinates);

    // Set the radius as the distance between the midpoint and the furthest point
    for(auto *v : data) {
        distance = get_distance(v->get_datum(), pivot);
        if(distance > radius) radius = distance;
    }

    if(number_of_points <= leaf_size) {
        child_1 = this;
        child_2 = this;
        contained_points = data;
    } else { // Make additional balls to break up data further
        contained_points = {};
        std::vector<Vertex *> v_y {y}, v_z {z};
        for(auto *v : data) {
            if(v != y && v != z) {
                if(get_distance(v, y) < get_distance(v, z)) {
                    v_y.push_back(v);
                } else {
                    v_z.push_back(v);
                }
            }
        }
        child_1 = new Ball(v_y, distance_function, leaf_size);
        child_2 = new Ball(v_z, distance_function, leaf_size);
    }
}

void Ball::initialize_mrd_ball(std::vector<Vertex *> data, TString metric, int leaf_size, Ball *query_tree, int param_k) {
    is_mrd_tree = true;
    contained_component = 0;
    distance_to_component = std::numeric_limits<double>::infinity();
    distance_function = metric;
    number_of_points = static_cast<int>(data.size());
    radius = 0;

    Vertex *x {data.back()}, *y {}, *z {};
    double max_distance {}, distance {};

    // y is the furthest point in the ball from x
    for(auto *v : data) {
        if(!y) {
            y = v;
            max_distance = get_distance(v, x);
        } else {
            distance = get_distance(v, x);
            if(distance > max_distance) {
                y = v;
                max_distance = distance;
            }
        }
    }

    // z is the furthest point in the ball from y
    max_distance = 0;
    for(auto *v : data) {
        if(v != y) {
            if(!z) {
                z = v;
                max_distance = get_distance(z, y);
            } else {
                distance = get_distance(z, y);
                if(distance > max_distance) {
                    z = v;
                    max_distance = distance;
                }
            }
        }
    }

    // Define a midpoint
    std::vector<double> midpoint_coordinates {};
    if(!z) z = y; // Should only happen with one point in ball
    for(int i = 0; i < y->get_datum()->get_number_dimensions(); ++i) {
        midpoint_coordinates.push_back((y->get_datum()->get_coordinate(i) + z->get_datum()->get_coordinate(i)) / 2.0);
    }
    pivot = new Point(midpoint_coordinates);
    std::vector<Point *> nn {};
    query_tree->find_knn(pivot, param_k, nn);
    pivot->set_core_distance(get_standard_distance(pivot, nn.back()));
    // Set the radius as the distance between the midpoint and furthest point
    for(auto *v : data) {
        distance = get_distance(v->get_datum(), pivot);
        if(distance > radius) radius = distance;
    }

    if(number_of_points <= leaf_size) {
        child_1 = this;
        child_2 = this;
        contained_points = data;
    } else { // Make additional balls to further break up data
        contained_points = {};
        std::sort(data.begin(), data.end(), [this, y](Vertex *a, Vertex *b) { return this->get_distance(a, y) < this->get_distance(b, y);});
        int halfway_point {static_cast<int>(data.size()) / 2};
        std::vector<Vertex *> v_1 {}, v_2 {};

        for(int i = 0; i < static_cast<int>(data.size()); ++i) {
            if(i < halfway_point) {
                v_1.push_back(data.at(i));
            } else {
                v_2.push_back(data.at(i));
            }
        }

        child_1 = new Ball(v_1, metric, leaf_size, query_tree, param_k);
        child_2 = new Ball(v_2, metric, leaf_size, query_tree, param_k);
    }
}

void Ball::find_contained_points(std::vector<Vertex *> &v) {
    if(is_leaf()) {
        v.insert(v.end(), contained_points.begin(), contained_points.end());
    } else {
        get_child_1()->find_contained_points(v);
        get_child_2()->find_contained_points(v);
    }
}

Ball::Ball(std::vector<Point *> data, TString metric, int leaf_size) {
    std::vector<Vertex *> contained_vertices {};
    for(Point *p : data) contained_vertices.push_back(new Vertex(p));
    initialize_standard_ball(contained_vertices, metric, leaf_size);
}

Ball::Ball(std::vector<Vertex *> data, TString metric, int leaf_size) {
    initialize_standard_ball(data, metric, leaf_size);
}

Ball::Ball(std::vector<Point *> data, TString metric, int leaf_size, Ball *query_tree, int param_k) {
    std::vector<Vertex *> vertices {};
    for(auto *p : data) {
        vertices.push_back(new Vertex(p));
        vertices.back()->set_datum(p);
    }
    initialize_mrd_ball(vertices, metric, leaf_size, query_tree, param_k);
}

Ball::Ball(std::vector<Vertex *> data, TString metric, int leaf_size, Ball *query_tree, int param_k) {
    initialize_mrd_ball(data, metric, leaf_size, query_tree, param_k);
}

std::vector<Vertex *> Ball::get_contained_data() {
    std::vector<Vertex *> vertices {};
    find_contained_points(vertices);
    return vertices;
}

double Ball::get_distance(Point *a, Point *b) {
    if(is_mrd_tree) {
        return get_mutual_reachability_distance(a, b);
    } else {
        return get_standard_distance(a, b);
    }
}

double Ball::get_distance(Vertex *a, Vertex *b) {
    if(is_mrd_tree) {
        return get_mutual_reachability_distance(a->get_datum(), b->get_datum());
    } else {
        return get_standard_distance(a->get_datum(), b->get_datum());
    }
}

void Ball::set_metric(TString metric) {
    if(metric == "E" || metric == "Euclidean") {
        distance_function = "Euclidean";
    } else {
        std::cout << "Ball::set_metric() : metric not found, using Euclidean" << std::endl;
        distance_function = "Euclidean";
    }
}

void Ball::set_contained_components() {
    if(contained_component) {
        contained_component = contained_component->get_component();
        return;
    }

    if(!is_leaf()) {
        Ball *c0 {get_child_1()}, *c1 {get_child_2()};
        Vertex *v0 {c0->get_contained_component()}, *v1 {c1->get_contained_component()};
        c0->set_contained_components();
        c1->set_contained_components();

        if(v0 && v1) {
            if(v0 == v1) { // Both children are set to the same cluster
                contained_component = v0;
            } else {
                contained_component = 0;
            }
        } else {
            contained_component = 0;
        }
    } else {
        for(int i = 0; i < static_cast<int>(contained_points.size()); ++i) {
            for(int j = i + 1; j < static_cast<int>(contained_points.size()); ++j) {
                // Keep components from being assigned to this cluster
                if(contained_points.at(i)->get_component() != contained_points.at(j)->get_component()) return;
            }
        }
        contained_component = contained_points.back()->get_component();
        return;
    }
}

void Ball::find_knn(Point *target, int k, std::vector<Point *> &found_neighbors) {
    // found_neighbors is modified recursively

    // If there are already k found neighbors and everything in this ball is already
    // further than the furthest neighbor, ignore this ball with no modification to
    // the found_neighbors vector.
    if(static_cast<int>(found_neighbors.size()) >= k && 
      (get_distance(target, pivot) - radius > get_distance(target, found_neighbors.back()))) return;

    if(is_leaf()) {
        for(auto *v : contained_points) {
            double distance {get_distance(target, v->get_datum())};
            bool add_point {};
            if(static_cast<int>(found_neighbors.size()) < k) {
                add_point = (target != v->get_datum()); // Never add the target point to its list
            } else {
                add_point = (target != v->get_datum() && distance < get_distance(target, found_neighbors.back()));
            }

            if(add_point) {
                if(found_neighbors.empty()) {
                    found_neighbors.push_back(v->get_datum());
                } else {
                    found_neighbors.insert(
                        std::lower_bound(found_neighbors.begin(), found_neighbors.end(), v, 
                                         [this, target](Point *p1, Vertex *p2) { 
                                            return get_distance(p1, target) < get_distance(p1, target); 
                                         }), v->get_datum());
                }

                if(static_cast<int>(found_neighbors.size()) > k) found_neighbors.pop_back();
            } // if(add_point)
        } // for(auto *v : contained_points)
    } else {
        if(get_distance(target, child_1->pivot) > get_distance(target, child_2->pivot)) {
            std::swap(child_1, child_2);
        }
        child_1->find_knn(target, k, found_neighbors);
        child_2->find_knn(target, k, found_neighbors);
    }
    return;
}

void Ball::find_knn(Point *target, int k, std::vector<Vertex *> &found_neighbors) {
    // found_neighbors is modified recursively

    // If there are already k found neighbors and everything in this ball is already
    // further than the furthest neighbor, ignore this ball with no modification to
    // the found_neighbors vector.
    if(static_cast<int>(found_neighbors.size()) >= k && 
      (get_distance(target, pivot) - radius > get_distance(target, found_neighbors.back()->get_datum()))) return;

    if(is_leaf()) {
        for(auto *v : contained_points) {
            double distance {get_distance(target, v->get_datum())};
            bool add_point {};
            if(static_cast<int>(found_neighbors.size()) < k) {
                add_point = (target != v->get_datum()); // Never add the target point to its list
            } else {
                add_point = (target != v->get_datum() && distance < get_distance(target, found_neighbors.back()->get_datum()));
            }

            if(add_point) {
                if(found_neighbors.empty()) {
                    found_neighbors.push_back(v);
                } else {
                    found_neighbors.insert(
                        std::lower_bound(found_neighbors.begin(), found_neighbors.end(), v, 
                                         [this, target](Vertex *v1, Vertex *v2) { 
                                            return get_distance(v1->get_datum(), target) < get_distance(v2->get_datum(), target); 
                                         }), v);
                }

                if(static_cast<int>(found_neighbors.size()) > k) found_neighbors.pop_back();
            } // if(add_point)
        } // for(auto *v : contained_points)
    } else {
        if(get_distance(target, child_1->pivot) > get_distance(target, child_2->pivot)) {
            std::swap(child_1, child_2);
        }
        child_1->find_knn(target, k, found_neighbors);
        child_2->find_knn(target, k, found_neighbors);
    }
    return;
}

void Ball::find_radial_nn(Point *target, double search_radius, std::vector<Vertex *> &found_neighbors) {
    // found_neighbors is modified recursively

    // If everyhting in the ball is outside of the search radius, ignore the ball
    if((get_distance(target, pivot) - radius > search_radius)) return;

    if(is_leaf()) {
        for(auto *v : contained_points) {
            double distance {get_distance(target, v->get_datum())};
            bool add_point {(target != v->get_datum()) && (distance < search_radius)};
            if(add_point) {
                if(found_neighbors.empty()) {
                    found_neighbors.push_back(v);
                } else {
                    found_neighbors.insert(
                        std::lower_bound(found_neighbors.begin(), found_neighbors.end(), v, 
                                         [this, target](Vertex *v1, Vertex *v2) { 
                                            return get_distance(v1->get_datum(), target) < get_distance(v2->get_datum(), target);
                                         }), v);
                }
            }
        }
    } else {
        if(get_distance(target, child_1->pivot) > get_distance(target, child_2->pivot)) std::swap(child_1, child_2);
        child_1->find_radial_nn(target, search_radius, found_neighbors);
        child_2->find_radial_nn(target, search_radius, found_neighbors);
    }
    return;
}

void Ball::find_radial_nn(Point *target, double search_radius, std::vector<Point *> &found_neighbors) {
    // found_neighbors is modified recursively

    // If everything in the ball is outside of the search radius, ignore the ball
    if((get_distance(target, pivot) - radius > search_radius)) return;

    if(is_leaf()) {
        for(auto *v : contained_points) {
            double distance {get_distance(target, v->get_datum())};
            bool add_point {(target != v->get_datum()) && (distance < search_radius)};
            if(add_point) {
                if(found_neighbors.empty()) {
                    found_neighbors.push_back(v->get_datum());
                } else {
                    found_neighbors.insert(std::lower_bound(found_neighbors.begin(), found_neighbors.end(), v->get_datum(),
                                           [this, target](Point *p1, Point *p2) { 
                                                return get_distance(p1, target) < get_distance(p2, target); 
                                            }), v->get_datum());
                }
            }
        }
    } else {
        if(get_distance(target, child_1->pivot) > get_distance(target, child_2->pivot)) std::swap(child_1, child_2);
        child_1->find_radial_nn(target, search_radius, found_neighbors);
        child_2->find_radial_nn(target, search_radius, found_neighbors);
    }
    return;
}

void Ball::find_radial_unsorted_nn(Point *target, double search_radius, std::vector<Point *> &found_neighbors) {
    // found_neighbors is modified recursively

    // If everything in the ball is outside of the search radius, ignore the ball
    if(get_distance(target, pivot) - radius > search_radius) return;

    if(is_leaf()) {
        for(auto *v : contained_points) {
            double distance {get_distance(target, v->get_datum())};
            bool add_point {target != v->get_datum() && distance < search_radius};
            if(add_point) found_neighbors.push_back(v->get_datum());
        }
    } else {
        if(get_distance(target, child_1->pivot) > get_distance(target, child_2->pivot)) std::swap(child_1, child_2);
        child_1->find_radial_unsorted_nn(target, search_radius, found_neighbors);
        child_2->find_radial_unsorted_nn(target, search_radius, found_neighbors);
    }
    return;
}

void Ball::find_radial_unsorted_nn(Point *target, double search_radius, std::vector<Vertex *> &found_neighbors) {
    // found_neighbors is modified recursively

    // If everything in the ball is outside of the search radius, ignore the ball
    if(get_distance(target, pivot) - radius > search_radius) return;

    if(is_leaf()) {
        for(auto *v : contained_points) {
            double distance {get_distance(target, v->get_datum())};
            bool add_point {(target != v->get_datum()) && (distance < search_radius)};
            if(add_point) found_neighbors.push_back(v);
        }
    } else {
        if(get_distance(target, child_1->pivot) > get_distance(target, child_2->pivot)) std::swap(child_1, child_2);
        child_1->find_radial_unsorted_nn(target, search_radius, found_neighbors);
        child_2->find_radial_unsorted_nn(target, search_radius, found_neighbors);
    }
    return;
}

void Ball::find_foreign_nn(Vertex *target, Vertex *neighbor) {
    // neighbor is modified recursively

    double distance_to_ball {get_distance(target->get_datum(), pivot) - radius};
    if(neighbor) {
        if(get_distance(target, neighbor) < distance_to_ball) return;
    }
    
    if(is_leaf()) {
        for(auto *v : contained_points) {
            double distance {get_distance(target, v)};
            bool new_point {};
            if(!neighbor) {
                new_point = (target != v) && (target->get_component() != v->get_component());
            } else {
                new_point = (target != v) && (target->get_component() != v->get_component()) 
                            && (distance < get_distance(target, v));
            }

            if(new_point) neighbor = v;
        }
    } else {
        if(get_distance(target->get_datum(), child_1->pivot) > get_distance(target->get_datum(), child_2->pivot)) {
            std::swap(child_1, child_2);
        }
        child_1->find_foreign_nn(target, neighbor);
        child_2->find_foreign_nn(target, neighbor);
    }
    return;
}

void Ball::reset_distance_to_components() {
    set_distance_to_component(std::numeric_limits<double>::infinity());
    if(!(this->is_leaf())) {
        // Continue resetting distances
        this->get_child_1()->reset_distance_to_components();
        this->get_child_2()->reset_distance_to_components();
    }
}

//__________________________________________________________________________________________________________________________________

void MinimumSpanningTree::find_nearest_foreign_vertex(Vertex *target, Ball *query_ball, Vertex *result) {
    double distance_to_edge_of_ball {query_ball->get_distance(target->get_datum(), query_ball->get_pivot()) - query_ball->get_radius()};

    if(result) {
        double distance_to_result {query_ball->get_distance(target->get_datum(), result->get_datum())};
        if(distance_to_edge_of_ball > distance_to_result) return;
    }

    if(query_ball->is_leaf()) {
        // do nothing?
    }
}

bool MinimumSpanningTree::completely_connected(Ball *a, Ball *b) {
    Vertex *a_component {a->get_contained_component()}, *b_component {b->get_contained_component()};
    if(a_component && b_component && (a_component == b_component)) return true;
    return false;
}

void MinimumSpanningTree::add_edges(std::vector<Vertex *> &forest) {
    // For every vertex in the forest:
    //      find the component c of that vertex
    //      if adding the edge found by find_component_neighbors connects two components:
    //           add the edge to the edge set
    for(auto *c: forest) {
        Vertex *c_1 {c->get_component()};
        if(c_1->get_foreign_neighbor()) {
            Vertex *v_1 {c_1->get_merging_vertex()};
            Vertex *v_2 {c_1->get_foreign_neighbor()};
            Vertex *c_2 {v_2->get_component()};
            if(c_1 != c_2) { // v_1 and v_2 are in different components
                edge_set.push_back(new Edge(v_1->get_datum(), v_2->get_datum(), c_1->get_delta()));
                c_2->set_parent(c_1);
            }
            c_1->set_delta(std::numeric_limits<double>::infinity());
            c_1->set_foreign_neighbor(0);
            c_1->set_merging_vertex(0);
        }
    }
}

MinimumSpanningTree::MinimumSpanningTree(Ball *ball_tree) {
    // Implements the dual tree Boruvka algorithm to construct a MST
    edge_sorted = false;
    vertex_set = ball_tree->get_contained_data();
    std::cout << "MinimumSpanningTree::MinimumSpanningTree() : number of vertices - " << vertex_set.size() << std::endl;

    std::vector<Vertex *> forest {vertex_set};
    edge_set = {};
    int edge_set_size {static_cast<int>(vertex_set.size()) - 1}; // number of edges to add

    // On the first pass, every vertex is in its own component. Thus, the first pass of
    // find_component_neighbors is equivalent to doing nearest neighbor searches for every point.
    std::cout << "MinimumSpanningTree::MinimumSpanningTree() : first pass of finding foreign nearest neighbors" << std::endl;
    std::vector<Vertex *> neighbors {};
    for(auto *v : vertex_set) {
        neighbors.clear();
        ball_tree->find_knn(v->get_datum(), 1, neighbors);
        v->set_foreign_neighbor(neighbors.back());
        v->set_merging_vertex(v);
        v->set_delta(ball_tree->get_distance(v->get_foreign_neighbor(), v->get_merging_vertex()));
        v->set_parent(v);
    }

    std::cout << "MinimumSpanningTree::MinimumSpanningTree() : adding first pass edges" << std::endl;
    add_edges(forest);
    std::cout << "MinimumSpanningTree::MinimumSpanningTree() : finding and adding the rest of the edges" << std::endl;

    while(static_cast<int>(edge_set.size()) < edge_set_size) {
        ball_tree->reset_distance_to_components();
        ball_tree->set_contained_components();
        find_component_neighbors(ball_tree, ball_tree);
        forest.erase(std::remove_if(forest.begin(), forest.end(), [this](Vertex *v) { return v->get_component() != v; }), forest.end());
        add_edges(forest);
    }
}

MinimumSpanningTree::MinimumSpanningTree(std::vector<Vertex *> v, std::vector<Edge *> e) {
    vertex_set = v;
    edge_set = e;
    edge_sorted = false;
};

// Edge *MinimumSpanningTree::get_nth_longest_edge();

bool MinimumSpanningTree::contains_edge(Edge *e) {
    for(auto *edge : edge_set) {
        if(e == edge) return true;
    }
    return false;
}

void MinimumSpanningTree::sort_edges() {
    // Sort from shortest to longest
    if(!edge_sorted) {
        std::sort(edge_set.begin(), edge_set.end(), [=](Edge *e1, Edge *e2) { return e1->get_length() < e2->get_length(); });
        edge_sorted = true;
    }
}

void MinimumSpanningTree::find_component_neighbors(Ball *a, Ball *b) {
    if(completely_connected(a, b)) {
        return;
    } else if(a->get_distance(a->get_pivot(), b->get_pivot()) 
            - a->get_radius() - b->get_radius() > a->get_distance_to_component()) {
                return;
    } else if(a->is_leaf() && b->is_leaf()) {
        for(auto *v_a : a->get_contained_data()) {
            Vertex *a_component {v_a->get_component()};
            a->set_distance_to_component(0);
            for(auto *v_b : b->get_contained_data()) {
                Vertex *b_component {v_b->get_component()};
                if(a_component != b_component) {
                    double distance {a->get_distance(v_a, v_b)};
                    if(distance < a_component->get_delta()) {
                        a_component->set_delta(distance);
                        a_component->set_foreign_neighbor(v_b);
                        a_component->set_merging_vertex(v_a);
                    }
                }
            }

            if(a->get_distance_to_component() < a_component->get_delta()) {
                a->set_distance_to_component(a_component->get_delta());
            }
        }
    } else {
        find_component_neighbors(a->get_child_1(), b->get_child_1());
        find_component_neighbors(a->get_child_2(), b->get_child_2());
        find_component_neighbors(a->get_child_2(), b->get_child_1());
        find_component_neighbors(a->get_child_1(), b->get_child_2());
        a->set_distance_to_component(std::max(a->get_child_1()->get_distance_to_component(), a->get_child_2()->get_distance_to_component()));
    }
}

//__________________________________________________________________________________________________________________________________

void HDBSCAN::set_extrema_values() {
    int number_dimensions {data.back()->get_number_dimensions()};
    if(static_cast<int>(minimum_values.size()) == 0 && static_cast<int>(maximum_values.size()) == 0) {
        for(int i = 0; i < number_dimensions; ++i) {
            minimum_values.push_back(std::numeric_limits<double>::infinity());
            maximum_values.push_back(-1.0 * std::numeric_limits<double>::infinity());
        }

        for(auto *p : data) {
            for(int i = 0; i < number_dimensions; ++i) {
                if(p->get_coordinate(i) < minimum_values.at(i)) minimum_values.at(i) = p->get_coordinate(i);
                if(p->get_coordinate(i) > maximum_values.at(i)) maximum_values.at(i) = p->get_coordinate(i);
            }
        }
    }
}

void HDBSCAN::dispose_mst() {
    if(minimum_spanning_tree) {
        delete minimum_spanning_tree;
        minimum_spanning_tree = 0;
    }
}

void HDBSCAN::dispose_ball_trees() {
    if(standard_metric_ball_root) {
        delete standard_metric_ball_root;
        standard_metric_ball_root = 0;
    }

    if(mutual_reachability_ball_root) {
        delete mutual_reachability_ball_root;
        mutual_reachability_ball_root = 0;
    }
}

void HDBSCAN::dispose_cluster_node_tree() {
    if(root_cluster) {
        delete root_cluster;
        root_cluster = 0;
    }
}
// Runs a k-th nearest neighbor query and saves the Euclidean nearest neighbors in a vector. The distance
// to the k-th NN is saved as the core distance. This function takes use of packages from alglib.
void HDBSCAN::initialize_ball_trees() {
    // Uses k parameter. Construct the standard metric ball with the contained data, specified metric
    // (usually Euclidean), and the maximum leaf size.
    if(!standard_metric_ball_root) {
        std::cout << "HDBSCAN::initialize_ball_trees() : constructing standard metric ball tree" << std::endl;
        standard_metric_ball_root = new Ball(data, metric, leaf_size);
        std::cout << "HDBSCAN::initialize_ball_trees() : assigning core distances" << std::endl;
        std::vector<Point *> knn {};
        for(auto *p : data) {
            knn.clear();
            standard_metric_ball_root->find_knn(p, k, knn);
            p->set_core_distance(knn.back()->get_mutual_reachability_distance(p));
        }
    } else {
        std::cout << "HDBSCAN::initialize_ball_trees() : standard_metric_ball_root already exists; it will be used"
                     "Core distances will not be updated." << std::endl;
    }

    // Once the spatial indexing tree is set up, run k nearest neighbor queries to determine core distances for the
    // mutual reachability ball
    if(!mutual_reachability_ball_root) {
        std::cout << "HDBSCAN::initialize_ball_trees() : constructing mutual reachability ball tree" << std::endl;
        // Use the standard metric ball to define the mutual reachability ball tree
        mutual_reachability_ball_root = new Ball(data, metric, leaf_size, standard_metric_ball_root, k);
    } else {
        std::cout << "HDBSCAN::initialize_ball_trees() : mutual_reachability_ball_root already exists; it will be used." << std::endl;
    }
}

void HDBSCAN::construct_minimum_spanning_tree() {
    // Calls constructor of minimum_spanning_tree which implements dual tree Boruvka algorithm to find a MST with use of a spatial indexing tree.
    if(!minimum_spanning_tree) {
        TStopwatch *clock {new TStopwatch()};
        clock->Start();
        minimum_spanning_tree = new MinimumSpanningTree(mutual_reachability_ball_root);
        std::cout << "HDBSCAN::construct_minimum_spanning_tree() : time to construct MST: " << clock->RealTime() << " seconds" << std::endl;
    } else {
        std::cout << "HDBSCAN::construct_minimum_spanning_tree() : minimum_spanning_tree already exists. A new tree will not be constructed." << std::endl;
    }
}

// Sort edges from short to long.
// for each edge:
//     if the data points have been assigned leaf clusters l1, l2:
//         s1, s2 = sizes of the two clusters (s1 smaller WLOG)
//         c1, c2 = l1.find, l2.find
//         if s1 < minimum_cluster_size: # small things coming together
//             assign all of the data points from c1 to c2
//             delete c1
//         else: # two large clusters merging
//             create a new cluster which is the parent of each
//     if one of the data points has been assigned a leaf cluster:
//         assign the point that doesn't have a cluster to the predefined cluster
//     if neither has been assigned a leaf cluster:
//         create a new cluster and assign both points to it

// l.find:
//     if l.parent is l:
//         return l
//     else:
//         l.parent.find
//         return l.parent

void HDBSCAN::construct_condensed_cluster_hierarchy() {
    if(!mutual_reachability_ball_root) initialize_ball_trees();
    if(!root_cluster) {
        minimum_spanning_tree->sort_edges();
        for(auto *e : minimum_spanning_tree->get_edge_set()) {
            merge_condensed_clusters(e->get_point(0), e->get_point(1), e->get_length());
        }
        std::cout << "HDBSCAN::construct_condensed_cluster_hierarchy() : assigning root cluster node" << std::endl;
        root_cluster = data.back()->get_cluster();
    } else {
        std::cout << "HDBSCAN::construct_condensed_cluster_hierarchy() : hierarchy already exists; a new hierarchy will not be calculated" << std::endl;
    }
}

void HDBSCAN::merge_condensed_clusters(Point *p1, Point *p2, double distance) {
    Node *c1 {p1->get_cluster()}, *c2 {p2->get_cluster()};
    // c_i is zero if p_i has not been assigned a cluster yet
    if(c1 && c2) {
        int s1 {static_cast<int>(c1->get_size())}, s2 {static_cast<int>(c2->get_size())};
        if(s2 < s1) { // Make c1 the smaller cluster
            std::swap(c1, c2);
            std::swap(s1, s2);
            std::swap(p1, p2);
        }

        if(s1 < minimum_cluster_size) {
            // Adding noise to either more noise or a large cluster
            for(auto *p : c1->get_contained_points()) {
                p->set_leaf_cluster(c2);
                c2->add_member_point(p);
            }
            c2->set_size(s2 + s1);
        } else {
            // Adding two clusters together
            c1->set_birth_lamdba(1.0 / distance);
            c2->set_birth_lamdba(1.0 / distance);
            c1->set_parent(new Node(c1, c2, distance));
            c2->set_parent(c1->get_parent());
        }
    } else if(!(c1 || c2)) {
        // Neither point has a defined leaf node yet.
        p1->set_leaf_cluster(new Node(p1));
        p2->set_leaf_cluster(p1->get_leaf_cluster());
        p1->get_cluster()->add_member_point(p2);
        p1->get_cluster()->set_size(2);
    } else {
        // Exactly one of the data points has a defined leaf node
        if(c1 && !c2) {
            std::swap(p1, p2);
            std::swap(c1, c2);
        }

        if(c2 && !c1) {
            p1->set_leaf_cluster(c2);
            c2->add_member_point(p1);
        } else {
            std::cout << "HDBSCAN::merge_condensed_clusters() : something has gone wrong" << std::endl;
        }
    }

    if(p1->get_lambda() < 0) p1->set_lambda(1.0 / distance);
    if(p2->get_lambda() < 0) p2->set_lambda(1.0 / distance);
}

// This is to do the excess of mass cluster selection. The logic breaks down into three main components:
// 1) if the cluster is a leaf node, it is the most stable of it or its children (of which it has none), then
// return itself back up the recursive chain; 2) if cluster is more stable than the most stable selection of
// clusters it owns, then nothing below it matters and return itself; 3) if the most stable selection of clusters
// below cluster are more stable than cluster, then that whole selection gets propagated back up the chain

std::vector<Node *> HDBSCAN::optimize_clusters_by_eom(Node *cluster) {
    if(cluster->get_child(0) == cluster) {
        // At leaf node
        if(cluster->get_size() < minimum_cluster_size) return {};

        // Set the death lambda to be the lambda which causes there to be less than minimum_cluster_size points
        std::vector<Point *> leaf_points {cluster->get_contained_points()};
        std::sort(leaf_points.begin(), leaf_points.end(), [](Point *p1, Point *p2) { return p1->get_lambda() > p2->get_lambda(); });
        cluster->set_death_lambda(leaf_points.at(minimum_cluster_size - 1)->get_lambda());
        cluster->set_default_stability();
        return {cluster};
    } else {
        double children_stability {};
        std::vector<Node *> stable_clusters {};
        // Find the most optimal cluster set on branch zero, one; add together
        for(auto *c : optimize_clusters_by_eom(cluster->get_child(0))) stable_clusters.push_back(c);
        for(auto *c : optimize_clusters_by_eom(cluster->get_child(1))) stable_clusters.push_back(c);
        for(auto *c : stable_clusters) children_stability += c->get_stability();

        if(alpha * children_stability > cluster->get_stability()) {
            return stable_clusters;
        } else {
            return {cluster};
        }
    }
}

std::vector<Node *> HDBSCAN::get_leaf_clusters(Node *cluster) {
    if(cluster->get_child(0) == cluster) {
        // At leaf node
        if(cluster->get_size() < minimum_cluster_size) {
            return {};
        } else {
            return {cluster};
        }
    } else {
        std::vector<Node *> leaf_clusters {};
        for(auto *c : get_leaf_clusters(cluster->get_child(0))) leaf_clusters.push_back(c);
        for(auto *c : get_leaf_clusters(cluster->get_child(1))) leaf_clusters.push_back(c);
        return leaf_clusters;
    }
}

void HDBSCAN::assign_clusters() {
    if(static_cast<int>(final_clusters.size()) == 0) {
        if(cluster_selection_method == "EOM") {
            if(allow_single_cluster) {
                final_clusters = optimize_clusters_by_eom(root_cluster);
                return;
            } else {
                final_clusters = {};
                for(auto *c : optimize_clusters_by_eom(root_cluster->get_child(0))) final_clusters.push_back(c);
                for(auto *c : optimize_clusters_by_eom(root_cluster->get_child(1))) final_clusters.push_back(c);
            }
        } else if(cluster_selection_method == "Leaf") {
            // Not currently equipped
        }

        label_clusters();
    } else {
        std::cout << "HDBSCAN::assign_clusters() : final clusters already calculated. A new set will not be calculated." << std::endl;
    }
}

void HDBSCAN::label_clusters() {
    for(auto *p : data) p->set_cluster_id(-1);

    for(int i = 0; i < static_cast<int>(final_clusters.size()); ++i) {
        Node *c {final_clusters.at(i)};
        for(auto *p : c->get_contained_points()) p->set_cluster_id(i);
    }
}

HDBSCAN::HDBSCAN(std::vector<Point *> &data, TString parameter_file_name) {
    this->data = data;

    set_parameters(parameter_file_name);
}

void HDBSCAN::set_parameters() {
    k = 15;
    minimum_cluster_size = 15;
    alpha = 1.0;
    leaf_size = 50.0;
    set_metric("Euclidean");
    scales = {};
    allow_single_cluster = false;
    minimum_spanning_tree = 0;
    standard_metric_ball_root = 0;
    mutual_reachability_ball_root = 0;
    root_cluster = 0;
    leaf_clusters = {};
    final_clusters = {};
    scaled = false;
    set_cluster_selection_method("EOM");
    set_extrema_values();
}

void HDBSCAN::set_parameters(TString parameter_file_name) {
    set_parameters();

    if(parameter_file_name != "") {
        std::fstream parameter_file(parameter_file_name, std::ios_base::in);
        parameter_file >> k >> minimum_cluster_size >> alpha >> leaf_size;
        parameter_file.close();
    }
}

void HDBSCAN::set_metric(TString metric) {
    if(metric == "E" || metric == "Euclidean") {
        this->metric = "E";
    } else {
        std::cout << "HDBSCAN::set_metric() : metric not found, using Euclidean" << std::endl;
        this->metric = "E";
    }

    dispose_mst();
    dispose_ball_trees();
    dispose_cluster_node_tree();
    final_clusters.clear();
}

void HDBSCAN::set_k(int k) {
    this->k = k;
    dispose_mst();
    dispose_ball_trees();
    dispose_cluster_node_tree();
    final_clusters.clear();
}

void HDBSCAN::set_ball_tree_leaf_size(double leaf_size) {
    this->leaf_size = leaf_size;
    dispose_ball_trees();
    final_clusters.clear();
}

void HDBSCAN::set_minimum_cluster_size(int minimum_size) {
    minimum_cluster_size = minimum_size;
    dispose_cluster_node_tree();
    for(auto *p : data) p->set_leaf_cluster(0);
    final_clusters.clear();
}

void HDBSCAN::set_alpha(double alpha) {
    this->alpha = alpha;
    final_clusters.clear();
}

void HDBSCAN::set_cluster_selection_method(TString metric) {
    if(metric == "EOM" || metric == "excess of mass") {
        cluster_selection_method = "EOM";
    } else if(metric == "Leaf") {
        cluster_selection_method = "Leaf";
    } else {
        std::cout << "HDBSCAN::set_cluster_selection_method() : selection method not found, using EOM" << std::endl;
        cluster_selection_method = "EOM";
    }

    final_clusters.clear();
}

void HDBSCAN::set_scales(std::vector<double> factors) {
    scales = factors;
    unscale_data();
    dispose_mst();
    dispose_ball_trees();
    dispose_cluster_node_tree();
    final_clusters.clear();
    scale_data();
}

void HDBSCAN::scale_data() {
    if(!scaled) {
        if(static_cast<int>(data.size()) == 0) return;

        // If scales are not set, set them to be one
        if(static_cast<int>(scales.size()) == 0) {
            for(int i = 0; i < data.back()->get_number_dimensions(); ++i) scales.push_back(1.0);
        }

        for(auto *p : data) {
            std::vector<double> new_coordinates {};
            for(int i = 0; i < static_cast<int>(scales.size()); ++i) {
                if(scales.at(i) <= 0) {
                    new_coordinates.push_back(p->get_coordinate(i));
                } else if(maximum_values.at(i) != minimum_values.at(i)) {
                    new_coordinates.push_back((p->get_coordinate(i) - minimum_values.at(i)) / (maximum_values.at(i) - minimum_values.at(i)) * scales.at(i));
                } else {
                    new_coordinates.push_back(0.0);
                }
            }
            p->set_point(new_coordinates);
        }
        scaled = true;
    }
}

void HDBSCAN::unscale_data() {
    if(scaled) {
        for(auto *p : data) {
            std::vector<double> new_coordinates {};
            for(int i = 0; i < p->get_number_dimensions(); ++i) {
                if(scales.at(i) > 0) {
                    new_coordinates.push_back(p->get_coordinate(i) * (maximum_values.at(i) - minimum_values.at(i)) / scales.at(i) + minimum_values.at(i));
                } else {
                    new_coordinates.push_back(p->get_coordinate(i));
                }
            }
            p->set_point(new_coordinates);
        }
        scaled = false;
    }
}

void HDBSCAN::print_parameters() {
    std::cout << std::endl;
    std::cout << "----------HDBSCAN::print_parameters()----------" << std::endl;
    std::cout << "K: " << k << std::endl;
    std::cout << "Minimum Cluster Size: " << minimum_cluster_size << std::endl;
    std::cout << "Alpha: " << alpha << std::endl;
    std::cout << "Leaf Size: " << leaf_size << std::endl;
    std::cout << "Allow Single Cluster: " << std::boolalpha << allow_single_cluster << std::endl;
    std::cout << "Cluster Selection Method: " << cluster_selection_method << std::endl;
    std::cout << std::endl;

    std::cout << "Scales: ";
    for(auto s : scales) std::cout << s << " ";
    std::cout << std::endl << std::endl;

    std::cout << "Variable Names: " << std::endl;
    for(auto n : variable_names) std::cout << n << " ";
    std::cout << std::endl << std::endl;

    std::cout << "Metric: " << metric << std::endl;
    std::cout << "Standard Metric Ball Root: " << standard_metric_ball_root << std::endl;
    std::cout << "Mutual Reachability Ball Root: " << mutual_reachability_ball_root << std::endl;
    std::cout << "Minimum Spanning Tree: " << minimum_spanning_tree << std::endl;
    std::cout << "Root Cluster: " << root_cluster << std::endl;
    std::cout << "Scaled: " << scaled << std::endl;
    std::cout << std::endl;
}

std::vector<double> HDBSCAN::get_scaled_coordinates(Point *p) {
    std::vector<double> new_coordinates {};
    for(int i = 0; i < static_cast<int>(scales.size()); ++i) {
        if(scales.at(i) <= 0) {
            new_coordinates.push_back(p->get_coordinate(i));
        } else if(maximum_values.at(i) != minimum_values.at(i)) {
            new_coordinates.push_back((p->get_coordinate(i) - minimum_values.at(i)) / 
                                      (maximum_values.at(i) - minimum_values.at(i)) * scales.at(i));
        } else {
            new_coordinates.push_back(0.0);
        }
    }
    return new_coordinates;
}

void HDBSCAN::fit() {
    TStopwatch *clock {new TStopwatch()};
    std::cout << "HDBSCAN::fit() : performing fit on " << static_cast<int>(data.size()) << " points" << std::endl;
    std::cout << "HDBSCAN::fit() : scaling data" << std::endl;
    scale_data();
    std::cout << "HDBSCAN::fit() : finding nearest neighbors" << std::endl;
    initialize_ball_trees();
    std::cout << "HDBSCAN::fit() : constructing the minimum spanning tree" << std::endl;
    construct_minimum_spanning_tree();
    std::cout << "HDBSCAN::fit() : constructing the condensed cluster hierarchy" << std::endl;
    construct_condensed_cluster_hierarchy();
    std::cout << "HDBSCAN::fit() : assigning clusters" << std::endl;
    assign_clusters();
    std::cout << "HDBSCAN::fit() : labeling spatial points" << std::endl;
    label_clusters();
    std::cout << "HDBSCAN::fit() : calculated clusters in " << clock->RealTime() << " seconds" << std::endl;
    unscale_data();
}

std::vector<double> HDBSCAN::get_unscaled_coordinates(Point *p) {
    std::vector<double> unscaled_coordinates {};
    for(int i = 0; i < p->get_number_dimensions(); ++i) {
        if(scales.at(i) > 0) {
            unscaled_coordinates.push_back(p->get_coordinate(i) * (maximum_values.at(i) - minimum_values.at(i)) 
                                           / scales.at(i) - minimum_values.at(i));
        } else {
            unscaled_coordinates.push_back(p->get_coordinate(i));
        }
    }
    return unscaled_coordinates;
}

#endif