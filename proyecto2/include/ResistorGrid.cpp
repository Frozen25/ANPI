//
// Created by ubuntu on 01/10/18.
//

#include "ResistorGrid.hpp"

namespace anpi {
  
  ResistorGrid::ResistorGrid() {
    //TODO: implementar funcion
  }
  
  ResistorGrid::~ResistorGrid() {
    //TODO: implementar funcion
  }
  
  std::size_t ResistorGrid::nodesToIndex(const std::size_t row1,
                                         const std::size_t col1,
                                         const std::size_t row2,
                                         const std::size_t col2) {
    //TODO: implementar funcion
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
    anpi::Matrix<float,std::allocator<float> > amapTmp(map.rows,
                                                       map.cols,
                                                       map.ptr<float>());
    // And transform it to a SIMD-enabled matrix
    anpi::Matrix<float> amap(amapTmp);
    
    cv::waitKey();

    return true;
  }
  
  bool ResistorGrid::navigate(const indexPair &nodes) {
    //TODO: implementar funcion
  }
}