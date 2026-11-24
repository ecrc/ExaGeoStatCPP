#include <helpers/LocationPatterns.hpp>
#include <cassert>
#include <vector>
#include <iostream>
using exageostat::helpers::LocationPatterns;
int main(){
 for(auto p:{"regular","chessboard","left-bottom","clusters"}){
  std::vector<double>x(100000),y(100000),z(100000),w(100000);
  std::srand(42); LocationPatterns::Fill(100000,x.data(),y.data(),p);
  std::srand(42); LocationPatterns::Fill(100000,z.data(),w.data(),p);
  assert(x==z && y==w); double mean=0;
  for(int i=0;i<100000;i++) {assert(x[i]>=0 && x[i]<1 && y[i]>=0 && y[i]<1);mean+=x[i];
   if(std::string(p)=="chessboard")assert((int(8*x[i])+int(8*y[i]))%2==0);

  }
  mean/=100000;assert(std::abs(mean-(std::string(p)=="left-bottom"?1./3:.5))<.01);
 }
 bool threw=false;try{LocationPatterns::Validate("satellite");}catch(const std::invalid_argument&){threw=true;}assert(threw);
 for(auto pattern:{"regular","chessboard","left-bottom","clusters"}) {
  const int n=90000; std::vector<double>x(n),y(n),z(n),a(n),b(n),c(n);
  std::srand(42);LocationPatterns::Fill3D(n,x.data(),y.data(),z.data(),pattern);
  std::srand(42);LocationPatterns::Fill3D(n,a.data(),b.data(),c.data(),pattern);
  assert(x==a && y==b && z==c); double means[3]={0,0,0};
  for(int i=0;i<n;i++) {
   assert(x[i]>=0 && x[i]<1 && y[i]>=0 && y[i]<1 && z[i]>=0 && z[i]<1);
   if(std::string(pattern)=="chessboard")assert((int(8*x[i])+int(8*y[i])+int(8*z[i]))%2==0);
   means[0]+=x[i];means[1]+=y[i];means[2]+=z[i];
  }
  const double target=std::string(pattern)=="left-bottom"?1./3:.5;
  for(double sum:means)assert(std::abs(sum/n-target)<.015);
 }
 std::cout<<"All pattern checks passed\n";
}
