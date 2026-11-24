#ifndef EXAGEOSTAT_LOCATION_PATTERNS_HPP
#define EXAGEOSTAT_LOCATION_PATTERNS_HPP
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
namespace exageostat::helpers {
// Paper-inspired sampling rules; not a reconstruction of the competition scripts.
struct LocationPatterns {
    static std::string Validate(const std::string &s) {
        if(s=="jittered-grid" || s=="regular" || s=="chessboard" || s=="left-bottom" || s=="clusters") return s;
        throw std::invalid_argument("Unknown location pattern. Use jittered-grid, regular, chessboard, left-bottom, or clusters. Satellite locations require input data.");
    }
    static double Uniform() { return static_cast<double>(std::rand())/(static_cast<double>(RAND_MAX)+1.0); }
    template<class T> static void Fill(int n, T *x, T *y, const std::string &pattern) {
        Validate(pattern);
        if(n<0 || (n && (!x || !y))) throw std::invalid_argument("Invalid location storage.");
        if(pattern=="jittered-grid") throw std::invalid_argument("Use the existing jittered-grid generator.");
        const double centers[5][2]={{.2,.2},{.2,.8},{.5,.5},{.8,.2},{.8,.8}};
        for(int i=0;i<n;++i) {
            double a,b;
            if(pattern=="clusters") {
                int c=static_cast<int>(5*Uniform());
                // Uniform disks, radius .12, entirely inside the unit square.
                double radius=.12*std::sqrt(Uniform()), angle=6.283185307179586*Uniform();
                a=centers[c][0]+radius*std::cos(angle); b=centers[c][1]+radius*std::sin(angle);
            } else {
                do { a=Uniform(); b=Uniform(); }
                while(pattern=="chessboard" && (static_cast<int>(8*a)+static_cast<int>(8*b))%2!=0);
                if(pattern=="left-bottom") { a*=a; b*=b; }
            }
            x[i]=static_cast<T>(a); y[i]=static_cast<T>(b);
        }
    }
};
}
#endif
