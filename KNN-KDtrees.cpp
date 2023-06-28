#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>


using namespace std;

// Define a point in 3D space
struct Point3D {
    double x, y, z;
    Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

// Define a node in the KD tree
struct KDNode {
    Point3D point;
    KDNode* left;
    KDNode* right;
    KDNode(Point3D p) : point(p), left(nullptr), right(nullptr) {}
};

// Define a comparator for the priority queue
struct QueueComparator {
    bool operator() (const pair<double, KDNode*>& lhs, const pair<double, KDNode*>& rhs) const {
        return lhs.first < rhs.first;
    }
};

// Define a function to calculate the distance between two points in 3D space
double euclideanDistance(Point3D p1, Point3D p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return (dx*dx + dy*dy + dz*dz);
}

// Define a recursive function to build the KD tree
KDNode* buildKDTree(vector<Point3D>& points, int depth) {
    if (points.empty()) {
        return nullptr;
    }

    int axis = depth % 3;

    auto compare = [axis](Point3D p1, Point3D p2) {
        if (axis == 0) return p1.x < p2.x;
        if (axis == 1) return p1.y < p2.y;
        return p1.z < p2.z;
    };

    auto begin = points.begin();
    auto end = points.end();
    auto median = begin + (end - begin) / 2;

    nth_element(begin, median, end, compare);

    KDNode* node = new KDNode(*median);

    vector<Point3D> leftPoints(begin, median);
    node->left = buildKDTree(leftPoints, depth + 1);

    vector<Point3D> rightPoints(median+1, end);
    node->right = buildKDTree(rightPoints, depth + 1);

    return node;
}

// Define a function to search for the k nearest neighbors in the KD tree
void searchKDTree(KDNode* node, Point3D query, priority_queue<pair<double, KDNode*>, vector<pair<double, KDNode*>>, QueueComparator> &pq, int k, int i) {

    if (node==nullptr) return;

    double curr_dist = euclideanDistance(node->point, query);


    if (pq.size() < k) {
        pq.push({curr_dist, node});

    } else if (curr_dist < pq.top().first) {
        pq.pop();
        pq.push({curr_dist, node});
    }

    double perp_; 

    if (i==0) perp_ = node->point.x - query.x;
    else if (i==1) perp_ = node->point.y - query.y;
    else perp_ = node->point.z - query.z;

    i = (i+1) % 3;

    //cout<<perp_ << "and D=" << curr_dist << " jhlj " << pq.size() <<endl;

    if (pow(perp_,2)<pq.top().first){
        searchKDTree(node->left, query, pq, k, i);
        searchKDTree(node->right, query, pq, k, i);
    }

    else {
        if (perp_<0) {
            searchKDTree(node->right, query, pq, k, i);
            }
        else searchKDTree(node->left, query, pq, k, i);
    }

    return;
}

/*
int main() {

    priority_queue<pair<double, KDNode*>, vector<pair<double, KDNode*>>, QueueComparator> pq;

    // Create some sample points
    vector<Point3D> points = {Point3D(0,0,0)};

    // Add more sample points
    points.push_back(Point3D(2.0,3.0,1.0));
    points.push_back(Point3D(4.0, 7.0, 8.0));
    points.push_back(Point3D(8.0, 1.0, 5.0));
    points.push_back(Point3D(5.0, 4.0, 2.0));
    
    // Build the KD tree
    KDNode* root = buildKDTree(points, 0);
    
    // Define a query point
    Point3D query(3.0, 5.0, 2.0);
    
    // Search for the 3 nearest neighbors of the query point
    searchKDTree(root, query, pq, 2, 0);
    
    vector<Point3D> neighbors;
    while (!pq.empty()) {
        neighbors.push_back(pq.top().second->point);
        pq.pop();
    }
    reverse(neighbors.begin(), neighbors.end());
    
    // Print the results
    cout << "The " << neighbors.size() << " nearest neighbors of the query point are:" << endl;
    for (auto neighbor : neighbors) {
        cout << "(" << neighbor.x << ", " << neighbor.y << ", " << neighbor.z << ")" << endl;
    }
    

    cout << root->point.x << "," << root->point.y << "," << root->point.z << endl;
    cout << root->left->point.x << "," << root->left->point.y << "," << root->left->point.z << endl;
    cout << root->right->point.x << "," << root->right->point.y << "," << root->right->point.z << endl;

    cout << root->left->left->point.x << "," << root->left->left->point.y << "," << root->left->left->point.z << endl;
    //cout << root->left->right->point.x << "," << root->left->right->point.y << "," << root->left->right->point.z << endl;

    cout << root->right->left->point.x << "," << root->right->left->point.y << "," << root->right->left->point.z << endl;
    //cout << root->right->right->point.x << "," << root->right->right->point.y << "," << root->right->right->point.z << endl;

    return 0;
}
*/