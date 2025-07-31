# **DENOISING: Dynamic Enhancement and Noise Overcoming in Multimodal Neural Observations via High-density CMOS-based Biosensors**

<p align="center">
  <img width="1090" height="614" src="https://github.com/HayderAminLab/DENOISING/assets/158823360/1c5ad1e1-5746-44b4-8a2c-92ceecde429b">
</p>

 A computational framework designed to enhance the fidelity and clarity of large-scale bioelectrical signals by dynamically mitigating noise. It addresses the challenge of 
 pervasive noise in neural recordings from high-density microelectrode arrays, offering an adaptable solution for various noise types and signal dimensions. By improving the signal-to-noise ratio 
 and facilitating the separation of signal and noise components, DENOISING aids in the accurate analysis and interpretation 
 of complex neural connectivity and dynamics without relying on specific data transformations.

## **1. Contents**

  - README - current file
  - LICENSE - usage and redistribution conditions
  - DENOISING.py - main top-level function used to denoise multimodal datasets
  - LFP_Denoising.py - top-level function used for denoising detected-LFP events
  - Spikes_Denoising.py - top-level function used for denoising detected-spikes

## **2. Getting Started**

  - Download data from Zenodo and separate it into subfolders based on datatype i.e. LFP and Spike.
  - If using personal data, format folder structure as instructed for the provided data. 
  - Clone the repository locally.
       > git clone https://github.com/HayderAminLab/DENOISING.git
  - Go to DENOISING.py and work through the script as follows:
    1. Change srcfilepath to folder with dataset
    2. Indicate data type (Analysis_Item='LFPs') OR (Analysis_Item='Spikes')
    3. Change parameters based on dataset using LFP_Denoising.py OR Spikes_Denoising.py as needed.
    4. Use Raster_plot function to compare raw (denoising=False) and denoised (denoising=True) raster data.
    5. Refine denoised datsets through repetition of steps iii-iv if needed. 

## **3. Data**

The following datasets have been provided for using and test running the DENOISING script.
  - Samples of LFP-detected and Spike-detected Datasets with Noise (https://doi.org/10.5281/zenodo.13284452), made available under Creative Commons Attribution 4.0 International!

## **4. Requirements**

Python >= 3.7; all analyses and testing were performed using Python 3.7 within PyCharm V.2023.2

## **5. 📄 Citation & Associated Publication**

We kindly ask you to **📌 cite** our paper if you use our code in your research.

**📘 Publication:**  
Hu X, Emery BA, Khanzada S, Amin H. (2024).  
**DENOISING**: Dynamic enhancement and noise overcoming in multimodal neural observations via high-density CMOS-based biosensors.  
*Frontiers in Bioengineering and Biotechnology*, 12:1390108.  
[🔗 View online](https://www.frontiersin.org/articles/10.3389/fbioe.2024.1390108/full) • [📄 Download PDF](https://github.com/HayderAminLab/DENOISING/raw/main/Hu%20et%20al%202024.pdf)

<pre> @article{Hu2024DENOISING,
  author    = {Xin Hu, Brett Addison Emery, Shahrukh Khanzada, and Hayder Amin},
  title     = {{DENOISING}: Dynamic enhancement and noise overcoming in multimodal neural observations via high-density {CMOS}-based biosensors},
  journal   = {Frontiers in Bioengineering and Biotechnology},
  volume    = {12},
  pages     = {1390108},
  year      = {2024},
  doi       = {10.3389/fbioe.2024.1390108},
  url       = {https://www.frontiersin.org/articles/10.3389/fbioe.2024.1390108/full}
}
</pre>

## **6. 📬 Contact**

For questions about the 🧠 **`code`**, please [open an issue](https://github.com/HayderAminLab/DENOISING/issues) in this repository.

For questions about the 📄 **`paper`**, feel free to contact  
**✉️ [Dr.-Ing. Hayder Amin](mailto:hayder.amin@dzne.de)** 
