//
// Created by ubuntu on 01/10/18.
//

#include "ResistorGrid.hpp"

namespace anpi {
  
  ResistorGrid::ResistorGrid() {
    ResistorGrid::rawMap_(1,1,0);
    ResistorGrid::A_(1,1,0);
    ResistorGrid::b_(1,0);
  }
  
  std::size_t ResistorGrid::nodesToIndex(const std::size_t row1,
                                         const std::size_t col1,
                                         const std::size_t row2,
                                         const std::size_t col2,
                                         const std::size_t n) {
    if(col1 == col2){    //Vertical case

      size_t fnk, result;
      fnk = (n-1) + (n + 1)*row2 + n*row1;    //Where n is equal to newCol2;
      result = fnk - (n-col2);
      return result;

    }

    else{               //Horizontal case
      size_t fnk, result;
      fnk = ((col2+(n-row2))-1) + ((col2+(n-row2)) + 1)*row2 + (col2+(n-row2))*row1;   //Where (col2+(n-row2)) es equal to new Col2;
      result = fnk - (n-row2);
      return result;
    }
  }
  
  indexPair ResistorGrid::indexToNodes(const std::size_t idx) {
    //TODO: implementar funcion
  }
  
  bool ResistorGrid::build(const std::string filename) {
    // TODO: probar que las filas y las columnas no se cambien
    // Build the name of the image in the data path
    std::string mapPath = std::string( ANPI_DATA_PATH ) + filename;

    // Read the image using the OpenCV
    cv::Mat_<float> map;

    cv::imread(mapPath.c_str(),
               CV_LOAD_IMAGE_GRAYSCALE).convertTo(map,CV_32FC1);
    map /= 255.0f; // normalize image range to 0 .. 255

    // Convert the OpenCV matrix into an anpi matrix
    // We have to use the std::allocator to avoid an exact stride
    anpi::Matrix<float,std::allocator<float> > amapTmp(static_cast<const size_t>(map.rows),
                                                       static_cast<const size_t>(map.cols),
                                                       map.ptr<float>());
    // And transform it to a SIMD-enabled matrix
    anpi::Matrix<float> amap(amapTmp);
    
    ResistorGrid::rawMap_ = amap;

    return true;
  }
  
  bool ResistorGrid::navigate(const indexPair &nodes) {
    //TODO: implementar funcion
  }
}