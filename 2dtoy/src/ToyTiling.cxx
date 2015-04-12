#include "WireCell2dToy/ToyTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"

using namespace WireCell;

ToyTiling::ToyTiling()
{
}

ToyTiling::ToyTiling(WireCell::Slice slice,WireCellSst::GeomDataSource gds){
  WireCell::Channel::Group group = slice.group();
  for (int i=0;i!=group.size();i++){
    const WireCell::GeomWire *wire = gds.by_channel_segment(group.at(i).first,0);
    if (wire->plane() == static_cast<WireCell::WirePlaneType_t>(0)){
      wire_u.push_back(wire);
    }else if (wire->plane() == static_cast<WireCell::WirePlaneType_t>(1)){
      wire_v.push_back(wire);
    }else if (wire->plane() == static_cast<WireCell::WirePlaneType_t>(2)){
      wire_w.push_back(wire);
    }
  }

  for (int i=0;i!=wire_u.size();i++){
    for (int j=0;j!=wire_v.size();j++){
      
    }
  }
  
  //std::cout << wire_u.size() << " " << wire_v.size() << " " << wire_w.size() << std::endl;
  
}


ToyTiling::~ToyTiling()
{
  //delete all the cells
  for (int i=0;i!=cell_all.size();i++){
    cell_all[i] = 0;
  }
}

GeomWireSelection ToyTiling::wires(const GeomCell& cell) const
{
    return GeomWireSelection();
}
	
GeomCellSelection ToyTiling::cells(const GeomWire& wire) const
{
    return GeomCellSelection();
}

GeomCell* ToyTiling::cell(const GeomWireSelection& wires) const
{
    return 0;
}
