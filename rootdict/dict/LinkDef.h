#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedclasses;

#pragma link C++ namespace WireCell;

#pragma link C++ class WireCell::Interface-!;
#pragma link C++ class WireCell::IContext-!;
#pragma link C++ class WireCell::IWireDatabase-!;
#pragma link C++ class WireCell::IWireGeometry-!;
#pragma link C++ class WireCell::ICellTiling-!;


#pragma link C++ namespace units;
#pragma link C++ defined_in namespace units;
#pragma link C++ global units::*;

#pragma link C++ class WireCell::GeomCell;
#pragma link C++ class WireCell::MergeGeomCell;
#pragma link C++ class WireCell::MergeGeomCellSet;
#pragma link C++ class WireCell::GeomCellMap;
#pragma link C++ class WireCell::GeomCellSet;
#pragma link C++ class WireCell::CellChargeMap;
#pragma link C++ class WireCell::CellIndexMap;
// apparently, rootcling doesn't like sets
// #pragma link C++ class WireCell::GeomCellSet;
#pragma link C++ class WireCell::GeomCellSelection;

#pragma link C++ class WireCell::GeomWire;
#pragma link C++ class WireCell::GeomWireSet;
#pragma link C++ class WireCell::WireChargeMap;
#pragma link C++ class WireCell::WireIndexMap;
#pragma link C++ class WireCell::MergeGeomWire;
#pragma link C++ class WireCell::GeomWireMap;
#pragma link C++ class WireCell::GeomWireWireMap;
#pragma link C++ class WireCell::GeomWireWiresMap;
// apparently, rootcling doesn't like sets
// #pragma link C++ class WireCell::GeomWireSet;
#pragma link C++ class WireCell::GeomWireSelection;

#pragma link C++ class WireCell::MergeCellCluster;
#pragma link C++ class WireCell::GeomCluster;
#pragma link C++ class WireCell::GeomClusterSet;

#pragma link C++ class WireCell::ChargeSequence;
#pragma link C++ class WireCell::Trace;
#pragma link C++ class WireCell::Frame;
//#pragma link C++ class WireCell::WireCharge;
//#pragma link C++ class WireCell::WireChargeCollection;
#pragma link C++ class WireCell::Slice;

#pragma link C++ class WireCell::Point;
#pragma link C++ class WireCell::PointVector;

#pragma link C++ class WireCell::Vector;
#pragma link C++ class WireCell::VectorPair;

#pragma link C++ function WireCell::box_intersection;

#endif
