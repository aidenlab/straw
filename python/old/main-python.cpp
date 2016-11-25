/*
  The MIT License (MIT)
 
  Copyright (c) 2011-2016 Broad Institute, Aiden Lab
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <sstream>
#include "straw.h"
using namespace std;
using namespace boost::python;

struct Bins 
{
  list x;
  list y;
  list counts;
};

// Converts a C++ vector to a python list
template <class T>
list toPythonList(std::vector<T> vector) {
  typename std::vector<T>::iterator iter;
  list mylist;
  for (iter = vector.begin(); iter != vector.end(); ++iter) {
    mylist.append(*iter);
  }
  return mylist;
}

Bins straw_python(string argv)
{
  vector<int> x;
  vector<int> y;
  vector<float> counts;
  string norm, fname, chr1loc, chr2loc, unit, size;
  int binsize;
  istringstream iss(argv);
  iss >> norm;
  iss >> fname;
  iss >> chr1loc;
  iss >> chr2loc;
  iss >> unit;
  iss >> size;
  
  binsize=stoi(size);
  straw(norm, fname, binsize, chr1loc, chr2loc, unit, x, y, counts);
  Bins b;
  b.x = toPythonList(x);
  b.y = toPythonList(y);
  b.counts = toPythonList(counts);
  return b;
  //int length=x.size();
    //for (int i=0; i<length; i++) {
      //  printf("%d\t%d\t%.14g\n", x[i], y[i], counts[i]);   
    //}
}

BOOST_PYTHON_MODULE(straw_ext)
{
  using namespace boost::python;
  def("straw", straw_python);
  class_<Bins>("Bins")
    .def_readonly("x", &Bins::x)
    .def_readonly("y", &Bins::y)
    .def_readonly("counts", &Bins::counts);
}
