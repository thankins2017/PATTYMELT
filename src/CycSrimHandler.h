#ifndef CYCSRIMHANDLER__H
#define CYCSRIMHANDLER__H

#include "CycSrim.h"

// To add new elements for consideration
// /opt/new/include/cycapp/CycSrim.h

// Reads the form "CycSrim::EMaterial" string from file and returns the corresponding
// material. Unfortunately, new elements have to be added here as needed. It's not
// efficient or clean, but it's what we're working with.
CycSrim::EMaterials evaluate_cycsrim_material(const std::string &material) {
    if(material == "CycSrim::SrimMaterialCarbon") {
        return CycSrim::SrimMaterialCarbon;
    } else if(material == "CycSrim::SrimMaterialAluminum") {
        return CycSrim::SrimMaterialAluminum;
    } else if(material == "CycSrim::SrimMaterialCD2") {
        return CycSrim::SrimMaterialCD2;
    } else if(material == "CycSrim::SrimMaterialSn") {
        return CycSrim::SrimMaterialSn;
    } else if(material == "CycSrim::SrimMaterialCu") {
        return CycSrim::SrimMaterialCu;
    } else if(material == "CycSrim::SrimMaterialAu") {
        return CycSrim::SrimMaterialAu;
    } else {
        std::invalid_argument("CycSrimHandler::evaluate_cycsrim_material : invalid material " + material + ". Returning CycSrim::SrimMaterialAu");
        return CycSrim::SrimMaterialAu;
    }
}

#endif