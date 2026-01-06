#include "image_compression.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp> // Sobel

using namespace compression;

ImageCompression::ImageCompression(std::string path, int n_points, int mode) : img_src(), d(1, 1) {
    img_src = cv::imread(path);
    
    if(img_src.empty()) {
        std::cerr << "Erro ao abrir imagem!" << std::endl;
        exit(1);
    }

    d = triangulation::Delaunay(img_src.rows, img_src.cols);
    std::srand(std::time(nullptr));


    d.add_point(Point(0, 0));
    d.add_point(Point(0, img_src.cols - 1));
    d.add_point(Point(img_src.rows - 1, 0));
    d.add_point(Point(img_src.rows - 1, img_src.cols - 1));
    
    int count = 4; 

    // sampling
    if (mode == 1) {
        std::cout << "Calculando mapa de bordas" << std::endl;
        
        cv::Mat gray, grad_x, grad_y, grad_magnitude;
        cv::cvtColor(img_src, gray, cv::COLOR_BGR2GRAY);
        cv::Sobel(gray, grad_x, CV_32F, 1, 0, 3);
        cv::Sobel(gray, grad_y, CV_32F, 0, 1, 3);
        cv::magnitude(grad_x, grad_y, grad_magnitude);
        cv::normalize(grad_magnitude, grad_magnitude, 0.0, 1.0, cv::NORM_MINMAX);

        float base_prob = 0.02f; 

        while (count < n_points) {
            int r = std::rand() % img_src.rows;
            int c = std::rand() % img_src.cols;
            
            float importance = grad_magnitude.at<float>(r, c);
            float probability = importance + base_prob;
            float dice = (float)std::rand() / RAND_MAX;

            if (dice < probability) {
                d.add_point(Point(r, c));
                count++;
                if (count % 1000 == 0) std::cout << "." << std::flush;
            }
        }
    } 
    // aleatório
    else {
        std::cout << " Escolhendo pontos uniformemente..." << std::endl;
        
        while (count < n_points) {
            int r = std::rand() % img_src.rows;
            int c = std::rand() % img_src.cols;
            
            d.add_point(Point(r, c));
            count++;
            
            if (count % 1000 == 0) std::cout << "." << std::flush;
        }
    }

    std::cout << "\nTriangulacao concluida." << std::endl;
}


double triangle_double_area(Point p1, Point p2, Point p3) {
    return (p1.x() - p3.x()) * (p2.y() - p3.y()) - (p2.x() - p3.x()) * (p1.y() - p3.y());
}

std::tuple<float, float, float> barycentric_coords(Point p1, Point p2, Point p3, Point p) {
    double area_total = triangle_double_area(p1, p2, p3);
    
    if (std::abs(area_total) < 1e-9) return {0.33f, 0.33f, 0.33f};

    double area1 = triangle_double_area(p, p2, p3); 

    double w1 = area1 / area_total;
    
    double area2 = triangle_double_area(p1, p, p3);
    double w2 = area2 / area_total;
    
    double w3 = 1.0 - w1 - w2;
    
    return { (float)w1, (float)w2, (float)w3 };
}

std::pair<float, float> linear_interp_line(Point p1, Point p2, Point p) {
    double dist_total = std::sqrt(std::pow(p1.x()-p2.x(), 2) + std::pow(p1.y()-p2.y(), 2));
    double dist_p1 = std::sqrt(std::pow(p.x()-p1.x(), 2) + std::pow(p.y()-p1.y(), 2));
    
    if (dist_total < 1e-9) return {0.5f, 0.5f};
    
    float w2 = (float)(dist_p1 / dist_total); 
    float w1 = 1.0f - w2;
    return {w1, w2};
}

cv::Scalar ImageCompression::interpolate(Point p) {
    triangulation::Location l = d.locate(p);
    
    switch (l.t) {
        case triangulation::Location::Type::VERTEX: {
            int x = std::min((int)p.x(), img_src.rows - 1);
            int y = std::min((int)p.y(), img_src.cols - 1);
            cv::Vec3b c = img_src.at<cv::Vec3b>(x, y);
            return cv::Scalar(c[0], c[1], c[2]);
        }
        case triangulation::Location::Type::LINE: {
            auto [p1, p2] = d.points(l.halfedge);
            auto [w1, w2] = linear_interp_line(p1, p2, p);
            
            cv::Vec3b c1 = img_src.at<cv::Vec3b>((int)p1.x(), (int)p1.y());
            cv::Vec3b c2 = img_src.at<cv::Vec3b>((int)p2.x(), (int)p2.y());
            
            return cv::Scalar(w1*c1[0] + w2*c2[0],
                              w1*c1[1] + w2*c2[1],
                              w1*c1[2] + w2*c2[2]);
        }
        case triangulation::Location::Type::FACE: {
            auto [p1, p2, p3] = d.points(l.face);
            
            auto [w1, w2, w3] = barycentric_coords(p1, p2, p3, p);
            
            cv::Vec3b c1 = img_src.at<cv::Vec3b>((int)p1.x(), (int)p1.y());
            cv::Vec3b c2 = img_src.at<cv::Vec3b>((int)p2.x(), (int)p2.y());
            cv::Vec3b c3 = img_src.at<cv::Vec3b>((int)p3.x(), (int)p3.y());
            
            return cv::Scalar(w1*c1[0] + w2*c2[0] + w3*c3[0],
                              w1*c1[1] + w2*c2[1] + w3*c3[1],
                              w1*c1[2] + w2*c2[2] + w3*c3[2]);
        }
        default: 
            return cv::Scalar(0, 0, 0);
    }
}


void ImageCompression::run() {
    img_compression = cv::Mat(img_src.rows, img_src.cols, CV_8UC3, cv::Scalar(0, 0, 0));
    
    double sum_sq_diff = 0.0;
    long total_pixels = img_src.rows * img_src.cols;

    std::cout << "Rasterizando e calculando RMSE..." << std::endl;

    for (int i = 0; i < img_compression.rows; i++) {
        if(i % 50 == 0) std::cout << "Linha " << i << " / " << img_compression.rows << std::endl;

        for (int j = 0; j < img_compression.cols; j++) {
            Point p(i, j);
            cv::Scalar color = interpolate(p);
            
            cv::Vec3b &pixel_out = img_compression.at<cv::Vec3b>(i, j);
            pixel_out = cv::Vec3b((uchar)color[0], (uchar)color[1], (uchar)color[2]);

            cv::Vec3b pixel_in = img_src.at<cv::Vec3b>(i, j);
            
            double diffB = pixel_in[0] - pixel_out[0];
            double diffG = pixel_in[1] - pixel_out[1];
            double diffR = pixel_in[2] - pixel_out[2];
            
            sum_sq_diff += (diffB*diffB + diffG*diffG + diffR*diffR);
        }
    }


    double mse = sum_sq_diff / (total_pixels * 3.0);
    double rmse = std::sqrt(mse);

    std::cout << "\nRMSE: " << rmse << std::endl;

    std::string output_filename = "output.png"; 
    cv::imwrite(output_filename, img_compression);
    std::cout << "Imagem salva em: " << output_filename << std::endl;

    cv::imshow("Original", img_src);
    cv::imshow("Compression", img_compression);
    cv::waitKey(0);
    cv::destroyAllWindows();
}