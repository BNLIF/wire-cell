#include "WireCell2dToy/ToyTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"
#include <cmath>

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

  float dis_u[2],dis_v[2],dis_w[2],
    dis_puv[4],dis_puw[4],dis_pwv[4];
  int k_save;

  for (int i=0;i!=wire_u.size();i++){
    dis_u[0] = gds.wire_dist(*wire_u[i])-gds.pitch(static_cast<WireCell::WirePlaneType_t>(0))/2.;
    dis_u[1] = gds.wire_dist(*wire_u[i])+gds.pitch(static_cast<WireCell::WirePlaneType_t>(0))/2.;
    for (int j=0;j!=wire_v.size();j++){
      dis_v[0] = gds.wire_dist(*wire_v[j])-gds.pitch(static_cast<WireCell::WirePlaneType_t>(1))/2.;
      dis_v[1] = gds.wire_dist(*wire_v[j])+gds.pitch(static_cast<WireCell::WirePlaneType_t>(1))/2.;
      
      //four vertices around
      PointVector puv;
      PointVector puw;
      PointVector pwv;
      
      PointVector pcell;
      int ncell = 1;

      puv.push_back(gds.crossing_point(dis_u[0],dis_v[0],static_cast<WireCell::WirePlaneType_t>(0),static_cast<WireCell::WirePlaneType_t>(1)));
      puv.push_back(gds.crossing_point(dis_u[0],dis_v[1],static_cast<WireCell::WirePlaneType_t>(0),static_cast<WireCell::WirePlaneType_t>(1)));
      puv.push_back(gds.crossing_point(dis_u[1],dis_v[1],static_cast<WireCell::WirePlaneType_t>(0),static_cast<WireCell::WirePlaneType_t>(1)));
      puv.push_back(gds.crossing_point(dis_u[1],dis_v[0],static_cast<WireCell::WirePlaneType_t>(0),static_cast<WireCell::WirePlaneType_t>(1)));
      
      
      for (int k=0;k!=4;k++){
	dis_puv[k] = gds.wire_dist(puv[k],static_cast<WireCell::WirePlaneType_t>(2));
	//std::cout << dis_puv[k] << " ";
      }
      //std::cout << std::endl;

      int flag = 0;

      for (int k=0;k!=wire_w.size();k++){
  	dis_w[0] = gds.wire_dist(*wire_w[k])-gds.pitch(static_cast<WireCell::WirePlaneType_t>(2))/2.;
  	dis_w[1] = gds.wire_dist(*wire_w[k])+gds.pitch(static_cast<WireCell::WirePlaneType_t>(2))/2.;	
	
	//std::cout << dis_w[0] << " " << dis_w[1] << std::endl;

	for (int m = 0;m!=4;m++){
	  if (dis_puv[m] > dis_w[0] && dis_puv[m] < dis_w[1]){
	    flag = 1;
	    pcell.push_back(puv[m]);
	  }
	}
	if (flag==1) break;
      }

      // initialize uw and vw points
      if (flag==1){
	puw.push_back(gds.crossing_point(dis_u[0],dis_w[0],static_cast<WireCell::WirePlaneType_t>(0),static_cast<WireCell::WirePlaneType_t>(2)));
	puw.push_back(gds.crossing_point(dis_u[0],dis_w[1],static_cast<WireCell::WirePlaneType_t>(0),static_cast<WireCell::WirePlaneType_t>(2)));
	puw.push_back(gds.crossing_point(dis_u[1],dis_w[1],static_cast<WireCell::WirePlaneType_t>(0),static_cast<WireCell::WirePlaneType_t>(2)));
	puw.push_back(gds.crossing_point(dis_u[1],dis_w[0],static_cast<WireCell::WirePlaneType_t>(0),static_cast<WireCell::WirePlaneType_t>(2)));
	for (int k=0;k!=4;k++){
	  dis_puw[k] = gds.wire_dist(puw[k],static_cast<WireCell::WirePlaneType_t>(1));
	  if (dis_puw[k] > dis_v[0] && dis_puw[k] < dis_v[1]){
	    pcell.push_back(puw[k]);
	  }
	}

	pwv.push_back(gds.crossing_point(dis_v[0],dis_w[0],static_cast<WireCell::WirePlaneType_t>(1),static_cast<WireCell::WirePlaneType_t>(2)));
	pwv.push_back(gds.crossing_point(dis_v[0],dis_w[1],static_cast<WireCell::WirePlaneType_t>(1),static_cast<WireCell::WirePlaneType_t>(2)));
	pwv.push_back(gds.crossing_point(dis_v[1],dis_w[1],static_cast<WireCell::WirePlaneType_t>(1),static_cast<WireCell::WirePlaneType_t>(2)));
	pwv.push_back(gds.crossing_point(dis_v[1],dis_w[0],static_cast<WireCell::WirePlaneType_t>(1),static_cast<WireCell::WirePlaneType_t>(2)));
	for (int k=0;k!=4;k++){
	  dis_pwv[k] = gds.wire_dist(pwv[k],static_cast<WireCell::WirePlaneType_t>(0));
	  if (dis_pwv[k] > dis_u[0] && dis_pwv[k] < dis_u[1]){
	    pcell.push_back(pwv[k]);
	  }
	}
	
	//order all the points by phi angle
	GeomCell *cell = new GeomCell(ncell,pcell);
	
	// check order 
	//	pcell = cell->boundary();
	// std::cout << "Cell Count: " << pcell.size() << " " << cell->cross_section() << std::endl;
	// for (int k=0;k!=pcell.size();k++){
	//   std::cout << pcell[k].y << " " << pcell[k].z << " " << std::atan2(pcell[k].z - cell->center().z, pcell[k].y-cell->center().y) << std::endl;
	// }

	cell_all.push_back(cell);
	ncell++;
	
	

      }
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
