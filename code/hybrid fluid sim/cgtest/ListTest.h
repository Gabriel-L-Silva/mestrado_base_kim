#include "core/List.h"
#include "core/ObjectList.h"
#include <iostream>

struct X
{
  int a, b;
};

inline std::ostream&
operator <<(std::ostream& os, const X& x)
{
  return os << "a:" << x.a << " b:" << x.b;
}

struct Vertex;

using VertexAllocator = cg::BlockAllocator<Vertex, 32>;

struct Vertex: public cg::SharedObject,
  public cg::ObjectListNode<Vertex, VertexAllocator>
{
  float x, y, z;

  Vertex()
  {
    x = y = z = 0;
  }

  Vertex(float x, float y, float z):
    x{x}, y{y}, z{z}
  {

  }
};

inline std::ostream&
operator <<(std::ostream& os, const Vertex& v)
{
  return os << "x:" << v.x << " y:" << v.y << " z:" << v.z;
}

int
testList()
{
  cg::List<int> ilist;
  cg::List<X> xlist;
  cg::ObjectList<Vertex> vlist;

  for (auto i = 1; i < 10; ++i)
  {
    auto f = float(i);

    ilist.add(i);
    xlist.add({i, i << 1});
    vlist.insert(new Vertex{f, f + 1, f + 2});
  }
  std::cout << "Integer list\n";
  for (auto i : ilist)
    std::cout << i << '\n';
  std::cout << "\nStructured type list (X)\n";
  for (const auto& x : xlist)
    std::cout << x << '\n';
  std::cout << "\nObject list (Vertex)\n";
  for (auto& v : vlist)
    std::cout << *v << '\n';
  for (auto v = vlist.begin(); v != vlist.end(); ++v)
    std::cout << *v << '\n';
  getchar();
  return 0;
}
