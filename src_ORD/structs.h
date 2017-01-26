#ifndef H_STRUCTS_H
#define H_STRUCTS_H


struct par {
    double r[3];
    double q;
    double m;
    double p[3];
    double f[3];
    double peng;
};

struct foreng {
    double f[3];
    double peng;
};

/* treecode node */
struct tnode {
    int numpar, ibeg, iend;
    double x_min, y_min, z_min;
    double x_max, y_max, z_max;
    double x_mid, y_mid, z_mid;
    double radius, sqradius, aspect;
    int level, num_children, exist_ms;
    double ***ms;
    struct tnode *child[8];      //Child is ptr to array of ptrs to tnode children
	int local;
};

/* reduced treecode node for MPI sending */
struct tnodesend {
    int numpar, ibeg, iend;
    double x_min, y_min, z_min;
    double x_max, y_max, z_max;
    double x_mid, y_mid, z_mid;
    double radius, sqradius, aspect;
    int level, num_children, exist_ms;
    // double ***ms;
    // struct tnode *child[8];      //Child is ptr to array of ptrs to tnode children
    int local;
};


#endif
