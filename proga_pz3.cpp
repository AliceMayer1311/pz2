#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <sstream>
#include <Windows.h>

class PGMImage {
private:
    int width, height, maxVal;
    std::vector<std::vector<int>> imageData;

    int clamp(int value, int min, int max) {
        if (value < min) return min;
        if (value > max) return max;
        return value;
    }

    // Функция для пропуска комментариев
    void skipComments(std::ifstream& file) {
        char c;
        while (file >> std::ws && file.peek() == '#') {
            std::string comment;
            std::getline(file, comment);
        }
    }

public:
    PGMImage() : width(0), height(0), maxVal(255) {}

    bool load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла: " << filename << std::endl;
            return false;
        }

        // Чтение магического числа
        std::string magic;
        file >> magic;
        if (magic != "P2" && magic != "P5") {
            std::cerr << "Неверный формат PGM файла: " << filename << " (магическое число: " << magic << ")" << std::endl;
            return false;
        }

        std::cout << "Формат: " << magic << std::endl;

        // Пропуск комментариев и чтение размеров
        skipComments(file);
        file >> width;
        skipComments(file);
        file >> height;
        skipComments(file);
        file >> maxVal;

        std::cout << "Размеры: " << width << "x" << height << ", MaxVal: " << maxVal << std::endl;

        if (width <= 0 || height <= 0) {
            std::cerr << "Неверные размеры изображения: " << width << "x" << height << std::endl;
            return false;
        }

        // Выделение памяти для данных изображения
        imageData.resize(height, std::vector<int>(width));

        if (magic == "P2") {
            // Текстовый PGM
            for (int i = 0; i < height; ++i) {
                for (int j = 0; j < width; ++j) {
                    if (!(file >> imageData[i][j])) {
                        std::cerr << "Ошибка чтения данных изображения в позиции (" << i << ", " << j << ")" << std::endl;
                        return false;
                    }
                }
            }
        } else if (magic == "P5") {
            // Бинарный PGM
            file.get(); // Пропуск одного байта (обычно это перевод строки)
            for (int i = 0; i < height; ++i) {
                for (int j = 0; j < width; ++j) {
                    unsigned char pixel;
                    file.read(reinterpret_cast<char*>(&pixel), 1);
                    if (!file) {
                        std::cerr << "Ошибка чтения бинарных данных изображения" << std::endl;
                        return false;
                    }
                    imageData[i][j] = pixel;
                }
            }
        }

        file.close();
        std::cout << "Изображение загружено: " << filename << " (" << width << "x" << height << ")" << std::endl;
        return true;
    }

    bool save(const std::string& filename) {
        if (width <= 0 || height <= 0) {
            std::cerr << "Неверные размеры изображения для сохранения" << std::endl;
            return false;
        }

        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Ошибка создания файла: " << filename << std::endl;
            return false;
        }

        file << "P2\n";
        file << width << " " << height << "\n";
        file << maxVal << "\n";

        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                file << imageData[i][j] << " ";
            }
            file << "\n";
        }

        file.close();
        std::cout << "Изображение сохранено: " << filename << std::endl;
        return true;
    }

    void addSaltPepperNoise(double noiseLevel = 0.05) {
        if (width <= 0 || height <= 0) {
            std::cerr << "Изображение не загружено для добавления шума" << std::endl;
            return;
        }

        std::random_device rd;
        std::mt19937 gen(rd());

        int totalPixels = width * height;
        int noisePixels = static_cast<int>(totalPixels * noiseLevel);

        for (int i = 0; i < noisePixels; ++i) {
            int x = std::uniform_int_distribution<>(0, width - 1)(gen);
            int y = std::uniform_int_distribution<>(0, height - 1)(gen);

            if (std::uniform_real_distribution<>(0.0, 1.0)(gen) < 0.5) {
                imageData[y][x] = 0;
            } else {
                imageData[y][x] = maxVal;
            }
        }

        std::cout << "Добавлен шум 'соль и перец' уровень: " << noiseLevel << std::endl;
    }

    void addGaussianNoise(double mean = 0, double sigma = 25) {
        if (width <= 0 || height <= 0) {
            std::cerr << "Изображение не загружено для добавления шума" << std::endl;
            return;
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> dis(mean, sigma);

        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                double noise = dis(gen);
                int newValue = static_cast<int>(imageData[i][j] + noise);
                imageData[i][j] = clamp(newValue, 0, maxVal);
            }
        }

        std::cout << "Добавлен гауссовский шум (mean=" << mean << ", sigma=" << sigma << ")" << std::endl;
    }

    void medianFilter(int kernelSize = 3) {
        if (width <= 0 || height <= 0) {
            std::cerr << "Изображение не загружено для фильтрации" << std::endl;
            return;
        }

        if (kernelSize % 2 == 0) {
            std::cerr << "Размер ядра должен быть нечетным" << std::endl;
            return;
        }

        std::vector<std::vector<int>> result = imageData;
        int offset = kernelSize / 2;

        for (int y = offset; y < height - offset; ++y) {
            for (int x = offset; x < width - offset; ++x) {
                std::vector<int> neighbors;

                for (int ky = -offset; ky <= offset; ++ky) {
                    for (int kx = -offset; kx <= offset; ++kx) {
                        neighbors.push_back(imageData[y + ky][x + kx]);
                    }
                }

                std::sort(neighbors.begin(), neighbors.end());
                result[y][x] = neighbors[neighbors.size() / 2];
            }
        }

        imageData = result;
        std::cout << "Применен медианный фильтр с ядром " << kernelSize << "x" << kernelSize << std::endl;
    }

    void meanFilter(int kernelSize = 3) {
        if (width <= 0 || height <= 0) {
            std::cerr << "Изображение не загружено для фильтрации" << std::endl;
            return;
        }

        std::vector<std::vector<int>> result = imageData;
        int offset = kernelSize / 2;

        for (int y = offset; y < height - offset; ++y) {
            for (int x = offset; x < width - offset; ++x) {
                int sum = 0;
                int count = 0;

                for (int ky = -offset; ky <= offset; ++ky) {
                    for (int kx = -offset; kx <= offset; ++kx) {
                        sum += imageData[y + ky][x + kx];
                        count++;
                    }
                }

                result[y][x] = sum / count;
            }
        }

        imageData = result;
        std::cout << "Применен усредняющий фильтр с ядром " << kernelSize << "x" << kernelSize << std::endl;
    }

    void gaussianFilter(int kernelSize = 3, double sigma = 1.0) {
        if (width <= 0 || height <= 0) {
            std::cerr << "Изображение не загружено для фильтрации" << std::endl;
            return;
        }

        std::vector<std::vector<int>> result = imageData;
        int offset = kernelSize / 2;

        std::vector<std::vector<double>> kernel(kernelSize, std::vector<double>(kernelSize));
        double sum = 0.0;

        for (int i = -offset; i <= offset; ++i) {
            for (int j = -offset; j <= offset; ++j) {
                double value = std::exp(-(i*i + j*j) / (2 * sigma * sigma));
                kernel[i + offset][j + offset] = value;
                sum += value;
            }
        }

        for (int i = 0; i < kernelSize; ++i) {
            for (int j = 0; j < kernelSize; ++j) {
                kernel[i][j] /= sum;
            }
        }

        for (int y = offset; y < height - offset; ++y) {
            for (int x = offset; x < width - offset; ++x) {
                double filteredValue = 0.0;

                for (int ky = -offset; ky <= offset; ++ky) {
                    for (int kx = -offset; kx <= offset; ++kx) {
                        filteredValue += imageData[y + ky][x + kx] * kernel[ky + offset][kx + offset];
                    }
                }

                result[y][x] = clamp(static_cast<int>(filteredValue), 0, maxVal);
            }
        }

        imageData = result;
        std::cout << "Применен фильтр Гаусса с ядром " << kernelSize << "x" << kernelSize << " (sigma=" << sigma << ")" << std::endl;
    }

    const std::vector<std::vector<int>>& getImageData() const {
        return imageData;
    }

    int getWidth() const { return width; }
    int getHeight() const { return height; }
    int getMaxVal() const { return maxVal; }
};

void testFilterEffectiveness(const PGMImage& original, const PGMImage& filtered) {
    if (original.getWidth() != filtered.getWidth() || original.getHeight() != filtered.getHeight()) {
        std::cerr << "Размеры изображений не совпадают для тестирования" << std::endl;
        return;
    }

    int height = original.getHeight();
    int width = original.getWidth();
    double mse = 0.0;
    double psnr = 0.0;

    const auto& origData = original.getImageData();
    const auto& filtData = filtered.getImageData();

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            double error = origData[i][j] - filtData[i][j];
            mse += error * error;
        }
    }
    mse /= (width * height);

    if (mse > 0) {
        psnr = 20 * std::log10(255.0 / std::sqrt(mse));
    }

    std::cout << "MSE: " << mse << std::endl;
    std::cout << "PSNR: " << psnr << " dB" << std::endl;
}

bool fileExists(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}

// Функция для создания тестового изображения если файлы не найдены
void createTestImage(const std::string& filename, int width = 100, int height = 100) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка создания тестового файла: " << filename << std::endl;
        return;
    }

    file << "P2\n";
    file << width << " " << height << "\n";
    file << "255\n";

    // Создание простого тестового изображения (градиент)
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int value = (i + j) % 256;
            file << value << " ";
        }
        file << "\n";
    }

    file.close();
    std::cout << "Создано тестовое изображение: " << filename << std::endl;
}

void processImage(const std::string& inputFile, const std::string& outputPrefix) {
    PGMImage image;

    std::cout << "\n=== Обработка " << inputFile << " ===" << std::endl;

    if (!image.load(inputFile)) {
        std::cerr << "Не удалось загрузить изображение: " << inputFile << std::endl;
        return;
    }

    // Проверка что изображение загружено корректно
    if (image.getWidth() <= 0 || image.getHeight() <= 0) {
        std::cerr << "Изображение имеет недопустимые размеры" << std::endl;
        return;
    }

    // Тест 1: Медианный фильтр для шума "соль и перец"
    PGMImage test1;
    test1 = image; // Используем оператор присваивания
    test1.addSaltPepperNoise(0.1);
    test1.save(outputPrefix + "_noisy_saltpepper.pgm");

    auto originalData = image.getImageData(); // Сохраняем оригинальные данные

    auto start = std::chrono::high_resolution_clock::now();
    test1.medianFilter(3);
    auto end = std::chrono::high_resolution_clock::now();
    test1.save(outputPrefix + "_filtered_median3.pgm");

    std::cout << "Результаты медианного фильтра:" << std::endl;
    testFilterEffectiveness(image, test1);
    std::cout << "Время: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " мс" << std::endl;

    // Тест 2: Усредняющий фильтр для гауссовского шума
    PGMImage test2;
    test2 = image;
    test2.addGaussianNoise(0, 20);
    test2.save(outputPrefix + "_noisy_gaussian.pgm");

    start = std::chrono::high_resolution_clock::now();
    test2.meanFilter(3);
    end = std::chrono::high_resolution_clock::now();
    test2.save(outputPrefix + "_filtered_mean3.pgm");

    std::cout << "Результаты усредняющего фильтра:" << std::endl;
    testFilterEffectiveness(image, test2);
    std::cout << "Время: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " мс" << std::endl;

    // Тест 3: Фильтр Гаусса для комбинированного шума
    PGMImage test3;
    test3 = image;
    test3.addSaltPepperNoise(0.05);
    test3.addGaussianNoise(0, 15);
    test3.save(outputPrefix + "_noisy_combined.pgm");

    start = std::chrono::high_resolution_clock::now();
    test3.gaussianFilter(5, 1.5);
    end = std::chrono::high_resolution_clock::now();
    test3.save(outputPrefix + "_filtered_gaussian5.pgm");

    std::cout << "Результаты фильтра Гаусса:" << std::endl;
    testFilterEffectiveness(image, test3);
    std::cout << "Время: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " мс" << std::endl;
}

int main() {
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    std::vector<std::string> inputFiles = {
        "image1.pgm", "image2.pgm", "image3.pgm", "image4.pgm", "image5.pgm"
    };

    std::cout << "=== ПРОГРАММА ДЛЯ УДАЛЕНИЯ ШУМОВ С ИЗОБРАЖЕНИЙ ===" << std::endl;

    // Проверка существования файлов и создание тестовых если нужно
    bool allFilesExist = true;
    for (const auto& file : inputFiles) {
        if (!fileExists(file)) {
            std::cout << "Файл не найден: " << file << " - создаем тестовый..." << std::endl;
            createTestImage(file);
            allFilesExist = false;
        }
    }

    if (!allFilesExist) {
        std::cout << "Были созданы тестовые файлы. Запустите программу снова." << std::endl;
        return 0;
    }

    // Обработка каждого изображения
    for (size_t i = 0; i < inputFiles.size(); ++i) {
        std::string outputPrefix = "result_image" + std::to_string(i + 1);
        processImage(inputFiles[i], outputPrefix);
        std::cout << std::string(50, '-') << std::endl;
    }

    std::cout << "\n=== ОБРАБОТКА ВСЕХ ИЗОБРАЖЕНИЙ ЗАВЕРШЕНА ===" << std::endl;

    return 0;
}
