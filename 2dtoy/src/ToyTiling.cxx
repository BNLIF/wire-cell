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
    wire_all.push_back(wire);
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
      int k_save;

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
	if (flag==1) {
	  k_save = k;
	  break;
	}
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
	
	// fill cellmap
	GeomWireSelection wiresel;
	wiresel.push_back(wire_u[i]);
	wiresel.push_back(wire_v[j]);
	wiresel.push_back(wire_w[k_save]);	
	cellmap[cell]=wiresel;

	//fillwiremap

	if (wiremap.find(wire_u[i]) == wiremap.end()){
	  //not found
	  GeomCellSelection cellsel;
	  cellsel.push_back(cell);
	  wiremap[wire_u[i]]=cellsel;
	}else{
	  //found
	  wiremap[wire_u[i]].push_back(cell);
	}
	
	if (wiremap.find(wire_v[j]) == wiremap.end()){
	  //not found
	  GeomCellSelection cellsel;
	  cellsel.push_back(cell);
	  wiremap[wire_v[j]]=cellsel;
	}else{
	  //found
	  wiremap[wire_v[j]].push_back(cell);
	}

	if (wiremap.find(wire_w[k_save]) == wiremap.end()){
	  //not found
	  GeomCellSelection cellsel;
	  cellsel.push_back(cell);
	  wiremap[wire_w[k_save]]=cellsel;
	}else{
	  //found
	  wiremap[wire_w[k_save]].push_back(cell);
	}
	
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
  if (cellmap.find(&cell) == cellmap.end()){
    //not found 
    return GeomWireSelection();
  }else{
    //found
    return cellmap.find(&cell)->second;
  }
    
}
	
GeomCellSelection ToyTiling::cells(const GeomWire& wire) const
{
  if (wiremap.find(&wire) == wiremap.end()){
    return GeomCellSelection();
  }else{
    return wiremap.find(&wire)->second;
  }
}

GeomCell* ToyTiling::cell(const GeomWireSelection& wires) const
{
  if (wires.size()!=3) return 0;
  const GeomWire *wire1 = wires[0];
  const GeomWire *wire2 = wires[1];
  const GeomWire *wire3 = wires[2];

  if (wire1->plane() == wire2->plane() ||
      wire1->plane() == wire3->plane() || 
      wire2->plane() == wire3->plane()) return 0;

  GeomCellSelection cells1 = cells(*wire1);
  GeomCellSelection cells2 = cells(*wire2);
  GeomCellSelection cells3 = cells(*wire3);
  
  for (int i = 0; i!=cells1.size(); i++){
    const GeomCell *cell1 = cells1[i];
    for (int j =0; j!=cells2.size(); j++){
      const GeomCell *cell2 = cells2[j];
      if (*cell1==*cell2){
	for (int k=0;k!=cells3.size();k++){
	  const GeomCell *cell3 = cells3[k];
	  if (*cell1 == *cell3){
	    //there is a problem here, not sure what to do 
	    //return cell1;
	    return 0;
	  }
	}
      }
    }
  }

  return 0;

}
