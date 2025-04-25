# Monitoring Single Vesicle Acetylcholine Release with 1D-CNN

This repository contains code and resources for the AI model developed in the study:  
**"Monitoring Single Vesicle Acetylcholine Storage and Sub-quantal Release by Living Neurons and Organoids"**.

We developed a **One-Dimensional Convolutional Neural Network (1D-CNN)** to automatically detect vesicle release events from electrophysiological signal traces. This model enables fast and reliable identification of sub-quantal release patterns based on the waveform of time-resolved electrical signals.

---

## 📌 Project Overview

In this study, a signal detection system was built to monitor acetylcholine release from single vesicles in neurons and organoids. Each event is recorded as a time-series signal (x-axis: time, y-axis: current), where vesicle release appears as characteristic peaks.

To automate the identification of these peaks, we trained a 1D-CNN model using labeled electrophysiological recordings, allowing the detection of both full and sub-quantal release events with high accuracy.

---

## 📁 Repository Structure

- `runs/`: Training output
- `best_model.pth`: Saved best model weights
- `Release_1DCNN_Model.ipynb`: Jupyter notebook for training and evaluation
- `scaler.pkl`: Saved scaler for normalization
- `*.pdf`: Model validation figures or reports
- `README.md`: Project introduction

## 🧠 Model Architecture

We use a 1D Convolutional Neural Network optimized for detecting peaks in noisy electrophysiological signals. Key components include:

1. Multiple convolutional layers with ReLU activation
2. Batch normalization and dropout for generalization
3. Final dense layer for binary classification (event / no-event)
