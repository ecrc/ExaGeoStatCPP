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
    template<class T> static void Fill3D(int n, T *x, T *y, T *z, const std::string &pattern) {
        Validate(pattern);
        if(n<0 || (n && (!x || !y || !z))) throw std::invalid_argument("Invalid 3D location storage.");
        if(pattern=="jittered-grid") throw std::invalid_argument("Use the existing jittered-grid generator.");
        for(int i=0;i<n;++i) {
            double a,b,c;
            if(pattern=="clusters") {
                if(Uniform()<.15) { a=Uniform(); b=Uniform(); c=Uniform(); }
                else {
                    // Eight corner clusters and a central cluster, with overlapping tails.
                    const int k=static_cast<int>(9*Uniform());
                    const double cx=k==8?.5:((k&1)?.8:.2);
                    const double cy=k==8?.5:((k&2)?.8:.2);
                    const double cz=k==8?.5:((k&4)?.8:.2);
                    const double sigma=k==8?.18:.13+.01*(k%3);
                    do {
                        double r=std::sqrt(-2*std::log(1-Uniform()));
                        double angle=6.283185307179586*Uniform();
                        a=cx+sigma*r*std::cos(angle); b=cy+sigma*r*std::sin(angle);
                        r=std::sqrt(-2*std::log(1-Uniform()));
                        angle=6.283185307179586*Uniform();
                        c=cz+sigma*r*std::cos(angle);
                    } while(a<0 || a>=1 || b<0 || b>=1 || c<0 || c>=1);
                }
            } else {
                do { a=Uniform(); b=Uniform(); c=Uniform(); }
                while(pattern=="chessboard" && (static_cast<int>(8*a)+static_cast<int>(8*b)+static_cast<int>(8*c))%2!=0);
                if(pattern=="left-bottom") { a*=a; b*=b; c*=c; }
            }
            x[i]=static_cast<T>(a); y[i]=static_cast<T>(b); z[i]=static_cast<T>(c);
        }
    }
    template<class T> static void Fill(int n, T *x, T *y, const std::string &pattern) {
        Validate(pattern);
        if(n<0 || (n && (!x || !y))) throw std::invalid_argument("Invalid location storage.");
        if(pattern=="jittered-grid") throw std::invalid_argument("Use the existing jittered-grid generator.");
        const double centers[5][2]={{.2,.2},{.2,.8},{.5,.5},{.8,.2},{.8,.8}};
        for(int i=0;i<n;++i) {
            double a,b;
            if(pattern=="clusters") {
                // Diffuse clusters with overlapping tails and 15% uniform background.
                // Reject out-of-domain draws instead of clipping onto the edges.
                if (Uniform() < .15) {
                    a=Uniform(); b=Uniform();
                } else {
                    const int c=static_cast<int>(5*Uniform());
                    const double sx[5]={.11,.14,.16,.12,.10};
                    const double sy[5]={.14,.10,.12,.10,.14};
                    do {
                        const double radius=std::sqrt(-2*std::log(1-Uniform()));
                        const double angle=6.283185307179586*Uniform();
                        a=centers[c][0]+sx[c]*radius*std::cos(angle);
                        b=centers[c][1]+sy[c]*radius*std::sin(angle);
                    } while(a<0 || a>=1 || b<0 || b>=1);
                }
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
