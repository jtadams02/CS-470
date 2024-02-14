#include <iostream>
#include <vector>

using namespace std;

int main(){
    vector<int> testVector = {1,2,3,4,5};
    for(int i=0;i<testVector.size();i++){
        cout << testVector[i] << endl;
    }
}