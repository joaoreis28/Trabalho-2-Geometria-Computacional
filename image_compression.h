#include <opencv2/core/mat.hpp>
#include <opencv2/core/types.hpp>
#include <string>

#include "triangulation.h"

namespace compression {

class ImageCompression {
   public:
    ImageCompression(std::string path, int n_points, int mode);
    void run();

   private:
    cv::Scalar interpolate(Point p);
    triangulation::Delaunay d;
    cv::Mat img_src;
    cv::Mat img_compression;
};

}  // namespace compression
