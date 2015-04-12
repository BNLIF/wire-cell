#ifndef WIRECELL_TOYTILING_H
#define WIRECELL_TOYTILING_H

#include "WireCellTiling/TilingBase.h"
#include "WireCellData/GeomWireCellMap.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/GeomDataSource.h"

namespace WireCell {
    /**
     *  A bogus tiling class that doesn't do anything.
     */
    class ToyTiling : public TilingBase { 
    public:
      ToyTiling();
      ToyTiling(WireCell::Slice slice,WireCellSst::GeomDataSource gds);
      virtual ~ToyTiling();
      
	GeomWireSelection wires(const GeomCell& cell) const;
	GeomCellSelection cells(const GeomWire& wire) const;
	virtual GeomCell* cell(const GeomWireSelection& wires) const;

	GeomCellSelection get_allcell(){ return cell_all;}
	GeomWireSelection get_allwire(){ return wire_all;}

	GeomWireSelection wire_u;
	GeomWireSelection wire_v;
	GeomWireSelection wire_w;

	GeomCellSelection cell_all;
	GeomWireSelection wire_all;
	
	GeomCellMap cellmap;
	GeomWireMap wiremap;

    };
}
#endif
