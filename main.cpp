#include <opencv2/core/types.hpp>
#include <opencv2/highgui.hpp>

#include "image_compression.h"

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::cout << "Usage: t2 /path/to/image.jpg" << std::endl;
        return 1;
    }

    std::string path = argv[1];
    int n_points = 0;
    int mode = 0;
    mode = std::stoi(argv[3]);

    n_points = std::stoi(argv[2]);
    compression::ImageCompression c(path, n_points, mode);

    c.run();
}
