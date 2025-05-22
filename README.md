
# Time Series Forecasting Project

This project focuses on modeling and forecasting a time series representing the *Industrial Production Index (IPI)*. The analysis involves several statistical models, with a particular emphasis on ARIMA modeling. The primary goal is to understand the underlying structure of the IPI data and to make accurate short-term forecasts.

## 📁 Project Structure

```
.
├── main.R             # Main script: complete analysis and modeling (in English)
├── FrenchVersion.R    # Full French version of the main analysis
├── hello.R            # Scratchpad file used for intermediate testing and exploration
├── IPI_202501.csv     # Dataset: Monthly Industrial Production Index
├── report_time_series_TZIEMI_TANGOUO.pdf         # Final report: contains written answers to specific questions
└── README.md          # Project documentation
```

## 📊 Dataset

* **File**: `IPI_202501.csv`
* **Description**: Contains monthly observations of the Industrial Production Index (IPI).
* **Objective**: Analyze and forecast the IPI using time series models.

## 🔍 Project Objectives

1. **Data Exploration & Visualization**
   Understand the structure, seasonality, and trends in the IPI data.

2. **Modeling**
   Fit time series models, specifically:

   * MA(2)
   * ARMA(1,2)
   * ARIMA(1,1,2)

3. **Model Comparison**
   Use likelihood ratio tests to compare nested models (e.g., MA(2) vs. ARMA(1,2)).

4. **Residual Diagnostics**
   Check for autocorrelation and model validity through residual analysis.

5. **Forecasting**
   Generate forecasts and confidence intervals for future observations.

6. **Confidence Regions**
   Construct and visualize joint prediction regions for forecasts using ellipses.

7. **Interpretation**
   Provide conclusions on model performance and prediction quality.

## 📌 Notes

* Some responses to project questions are provided directly in the **report** and are not included in `main.R`.
* The **FrenchVersion.R** script contains a fully translated version (in French) of the main analysis.
* The **hello.R** file was used as a scratch space for testing code snippets throughout the development process.

## ⚙️ How to Run

1. Open `main.R` in RStudio or any R environment.
2. Make sure to load the necessary libraries: `forecast`, `tseries`, `ggplot2`, `ellipse`, etc.
3. Ensure the file `IPI_202501.csv` is in the working directory.
4. Run the script from top to bottom to reproduce all analysis steps, including model fitting and forecasting.

## 📄 Requirements

* R (version ≥ 4.0.0)
* Required packages:

  ```r
  install.packages(c("forecast", "ggplot2", "tseries", "ellipse"))
  ```

## 🧠 Author

This project was developed by two students, as part of a coursework in time series analysis. It demonstrates key modeling steps, statistical testing, and forecast evaluation for real-world economic data.

