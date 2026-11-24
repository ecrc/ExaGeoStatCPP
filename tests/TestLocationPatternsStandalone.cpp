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
 std::cout<<"All pattern checks passed\n";
}
