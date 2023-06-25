use std::env;
use std::path::Path;

use num::complex;
use std::f64::consts::PI;

use image::GrayImage;
use image::imageops;
use image::imageops::FilterType;

fn dft_by_pixel(u: f64, v: f64, source_img: &GrayImage) -> f64 {
    let img_x = source_img.dimensions().0 as f64;
    let img_y = source_img.dimensions().1 as f64;
    let mut return_val = 0 as f64;
    for i in 0..(img_x as i32) {
        for j in 0..(img_y as i32) {
            let exponent = (u*f64::from(i) / img_x) + (v*f64::from(j) / img_y);
            let complex_part = complex::Complex::new(0.0, -2.0*PI*exponent);
            let image::Luma(source_val) = source_img.get_pixel(i as u32, j as u32);
            return_val += (source_val[0] as f64) * complex_part.exp().re;
        }
    }

    return return_val;
}


fn main() {
    let args: Vec<String> = env::args().collect();
    let img_path = Path::new(&args[1]);
    let img_filename = img_path.file_name().unwrap();
    let img_filestem_str = img_path.file_stem().unwrap().to_str().unwrap();
    let _img_extension_str = img_path.extension().unwrap().to_str().unwrap();
    
    let img = image::open(img_path).unwrap().into_luma8();
    let resize_img = imageops::resize(&img, 128, 128, FilterType::Gaussian);
    let img_x = resize_img.dimensions().0;
    let img_y = resize_img.dimensions().1;

    println!("Computing DFT of {}", img_filename.to_str().unwrap()); 
    let mut dft_values = Vec::new();
    let total_pixels = img_x * img_y;
    
    let mut current_pixels = 0;
    let mut max_dft = 0.0 as f64; 
    for u in 0..img_x {
        let mut row_dft_values = Vec::new();
        for v in 0..img_y {
            let dft = dft_by_pixel(u as f64, v as f64, &resize_img);
            row_dft_values.push(dft);
            
            max_dft = max_dft.max(dft);

            current_pixels += 1;
            print!("{} / {}, max = {:5.1}\r", current_pixels, total_pixels, max_dft);
        }
        dft_values.push(row_dft_values);
    }
    println!();

    
    println!("Logarithmic Scaling.");
    let mut fourier_img: GrayImage = image::ImageBuffer::new(img_x, img_y);
    let c = 255.0 / (1.0 + max_dft).ln();

    let mut current_pixels = 0;
    let total_pixels = img_x * img_y;
    for u in 0..img_x {
        for v in 0..img_y {
            let dft_value = dft_values[u as usize][v as usize];
            let log_dft = c * (dft_value.abs() + 1.0).ln();
            let pixel_val = num::clamp(log_dft, 0.0, 255.0) as u8;

            let pixel = fourier_img.get_pixel_mut(u, v);
            *pixel = image::Luma([pixel_val]);

            current_pixels += 1;
            print!("{} / {}\r", current_pixels, total_pixels);
        }
    }
    println!();


    println!("Swapping quadrants.");
    let img_x_half = img_x / 2;
    let img_y_half = img_y / 2;
    
    let mut current_pixels = 0;
    let total_pixels = img_x_half * img_y_half;
    for u in 0..img_x_half {
        for v in 0..img_y_half {
            // Swap a quadrant and c quandrant
            let a_pixel = image::Luma([fourier_img.get_pixel(u, v)[0]]);
            let c_pixel = image::Luma([fourier_img.get_pixel(u+img_x_half, v+img_x_half)[0]]);
            let mut mut_pixel = fourier_img.get_pixel_mut(u, v);
            *mut_pixel = c_pixel;
            mut_pixel = fourier_img.get_pixel_mut(u+img_x_half, v+img_x_half);
            *mut_pixel = a_pixel;

            // Swap b quandrant and d quandrant
            let b_pixel = image::Luma([fourier_img.get_pixel(u, v+img_y_half)[0]]);
            let d_pixel = image::Luma([fourier_img.get_pixel(u+img_x_half, v)[0]]);
            let mut mut_pixel = fourier_img.get_pixel_mut(u, v+img_y_half);
            *mut_pixel = d_pixel;
            mut_pixel = fourier_img.get_pixel_mut(u+img_x_half, v);
            *mut_pixel = b_pixel;

            current_pixels += 1;
            print!("{} / {}\r", current_pixels, total_pixels);
        }
    }
    println!();


    let fourier_img_stem = img_filestem_str.to_owned() + "_fourier";
    let fourier_img_path_str = fourier_img_stem + "." + "png";
    fourier_img.save(fourier_img_path_str).unwrap();
}
