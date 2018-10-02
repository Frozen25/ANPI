//
// Created by ubuntu on 01/10/18.
//

#include "ResistorGrid.hpp"

namespace anpi {
  std::size_t ResistorGrid::nodesToIndex(const std::size_t row1,
                           const std::size_t col1,
                           const std::size_t row2,
                           const std::size_t col2) {

  }

  // TODO: provar que las filas y las columnas no se cambien
  bool ResistorGrid::build(const std::string filename) {
    // Build the name of the image in the data path
    std::string mapPath = std::string( ANPI_DATA_PATH ) + filename;

    // Read the image using the OpenCV
    cv::Mat_<float> map;

    cv::imread(mapPath.c_str(),
               CV_LOAD_IMAGE_GRAYSCALE).convertTo(map,CV_32FC1);
    map /= 255.0f; // normalize image range to 0 .. 255

    // And create a window to show the image
    cv::namedWindow(mapPath,CV_WINDOW_NORMAL | CV_GUI_EXPANDED);
    cv::imshow(mapPath,map);

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
}