import java.lang.*;
import java.util.*;

public class SquareHough {
    private final static int gaussianWinSize = 7;
    private final static int gaussianWinMid = gaussianWinSize / 2;

    private final static int momentsWinSize = 7;
    private final static int momentsWinMid = momentsWinSize / 2;

    private final static int peaksWinSize = 19;
    private final static int peaksWinMid = peaksWinSize / 2;

    private static double squareSizeHalf;

    private static double deltaThetaDeg;
    private static double linePeakThreshold;
    private static double squarePeakThreshold;
    private static double backProjThreshold;

    private static boolean preapplySobel = false;
    private static boolean useCanny = false;
    private static String imgEdgesFileName = "DoG.pgm";

    private final static double[] sinDeg = new double[360];
    private final static double[] cosDeg = new double[360];

    public static void main(String[] args) {
        if (args.length != 7) {
            System.out.println("ERROR: number of arguments should be 7.");
            System.out.println("arg 1 - filename of the input PGM inputImage");
            System.out.println("arg 2 - size (length of sides) of the squares to be detected");
            System.out.println("arg 3 - delta theta");
            System.out.println("arg 4 - line accumulation threshold");
            System.out.println("arg 5 - square accumulation threshold");
            System.out.println("arg 6 - back projection threshold");
            System.out.println("arg 7 - L or E to indicate lines [DoG] or DoG of edges [SobelDoG] to be used");
            System.out.println("      - C will use a canny edge detector instead (extension)");
            return;
        }

        final Image inputImage = new Image();

        inputImage.ReadPGM(args[0]);

        if (inputImage.depth <= 0) {
            System.out.println("ERROR: image depth not valid, may not have found file");
            return;
        } else if (inputImage.width <= 0) {
            System.out.println("ERROR: image width not valid, may not have found file");
            return;
        } else if (inputImage.height <= 0) {
            System.out.println(("ERROR: image height not valid, may not have found file"));
            return;
        }

        double squareSize;

        try {
            squareSize = Double.parseDouble(args[1]);

            if (squareSize <= 0) {
                System.out.println("ERROR: arg 2 - size can't be smaller than or equal to 0");
                return;
            } else if (squareSize < 10) {
                System.out.println("WARN: arg 2 - chosen square size is small, may struggle for optimal results");
            }

            squareSizeHalf = squareSize / 2;
        } catch (NumberFormatException e) {
            System.out.println("ERROR: arg 2 not a real number");
            return;
        }

        try {
            deltaThetaDeg = Double.parseDouble(args[2]);

            if (deltaThetaDeg < 0) {
                System.out.println("ERROR: arg 3 - delta theta should be more than or equal to 0");
                return;
            } else if (squareSize < 2) {
                System.out.println("WARN: arg 3 - chosen delta theta is small, may struggle for optimal results");
            }
        } catch (NumberFormatException e) {
            System.out.println("ERROR: arg 3 not a real number");
            return;
        }

        try {
            linePeakThreshold = Double.parseDouble(args[3]);
        } catch (NumberFormatException e) {
            System.out.println("ERROR: arg 4 not a real number");
            return;
        }

        try {
            squarePeakThreshold = Double.parseDouble(args[4]);
        } catch (NumberFormatException e) {
            System.out.println("ERROR: arg 5 not a real number");
            return;
        }

        try {
            backProjThreshold = Double.parseDouble(args[5]);
        } catch (NumberFormatException e) {
            System.out.println("ERROR: arg 6 not a real number");
            return;
        }

        if (args[6].equals("C")) {
            useCanny = true;
            imgEdgesFileName = "Canny.pgm";
        } else if (args[6].equals("E")) {
            preapplySobel = true;
            imgEdgesFileName = "SobelDoG.pgm";
        } else if (!args[6].equals("L")) {
            System.out.println("ERROR: arg 7 - expected 'L' or 'E'");
            return;
        }

        for (int deg = 0; deg < 360; deg++) {
            sinDeg[deg] = Math.sin(Math.toRadians(deg));
            cosDeg[deg] = Math.cos(Math.toRadians(deg));
        }

        final Image imgEdges = imgEdgeDetection(inputImage);
        imgEdges.WritePGM(imgEdgesFileName);

        final double[][] lineAccumulator = houghTransformLine(imgEdges);
        final Image imgLineAccumulator = getImgLineAccumulator(lineAccumulator);
        imgLineAccumulator.WritePGM("accumulator.pgm");

        final ArrayList<Line> lines = new ArrayList<>();
        getPeaksLineAccumulator(lineAccumulator, lines);

        final ImagePPM imgLineOverlay = overlayLinesFromPeaks(inputImage, lines);
        imgLineOverlay.WritePPM("lines.ppm");

        final double[][][] squareAccumulator = houghTransformSquare(lineAccumulator);
        final ArrayList<Square> squares = new ArrayList<>();
        getPeaksSquareAccumulator(squareAccumulator, squares);

        final ArrayList<SquarePixels> squarePixelsList = getSquaresPixelsFromPeaks(inputImage, squares);
        final ImagePPM imgSquareOverlay = overlaySquaresFromPixels(inputImage, squarePixelsList);
        imgSquareOverlay.WritePPM("squares.ppm");

        final ImagePPM overlaySquareWithBackProjectionImage = overlaySquaresFromPixelsWithBackProj(inputImage, imgEdges, squarePixelsList);
        overlaySquareWithBackProjectionImage.WritePPM("backproject.ppm");
    }

    private static Image imgEdgeDetection(Image img) {
        final Image imgPreBlur;

        if (useCanny) {
            System.out.println("INFO: using canny edge detector");

            return cannyEdgeDetector(img);
        } else if (preapplySobel) {
            imgPreBlur = sobelEdgeDetector(img);

            rescaleImg(imgPreBlur);

            // Invert image
            for (int x = 0; x < imgPreBlur.width; x++) {
                for (int y = 0; y < imgPreBlur.height; y++) {
                    imgPreBlur.pixels[x][y] = imgPreBlur.depth / 2 - (imgPreBlur.pixels[x][y] - (imgPreBlur.depth / 2 + 1));
                }
            }

            System.out.println("INFO: applied inverted sobel edge detector");
        } else {
            imgPreBlur = new Image(img.depth, img.width, img.height);

            for (int x = 0; x < img.width; x++) {
                System.arraycopy(img.pixels[x], 0, imgPreBlur.pixels[x], 0, img.height);
            }
        }

        final Image imgSigma1Blur = gaussianBlur(imgPreBlur,1.0f);
        final Image imgSigma2Blur = gaussianBlur(imgPreBlur,2.0f);

        final Image difference = new Image(img.depth, img.width, img.height);

        for (int x = 0; x < img.width; x++) {
            for (int y = 0; y < img.height; y++) {
                difference.pixels[x][y] = Math.max(imgSigma2Blur.pixels[x][y] - imgSigma1Blur.pixels[x][y], 0);
            }
        }

        rescaleImg(difference);

        System.out.println("INFO: using difference of gaussian");

        return difference;
    }

    private static Image cannyEdgeDetector(Image img) {
        Image imgEdges = gaussianBlur(img, 1.2f);
        int[][] gradientsAngleDeg = new int[img.width][img.height];
        int[][] magnitude = new int[img.width][img.height];

        final int winSize = 3;
        final int winMid = winSize / 2;
        final int[][] kernel = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};

        int gradientsX;
        int gradientsY;

        for (int x = 0; x < img.width; x++) {
            for (int y = 0; y < img.height; y++) {
                gradientsX = 0;

                for (int xx = -winMid; xx < winMid + 1; xx++) {
                    for (int yy = -winMid; yy < winMid + 1; yy++) {
                        final int clampedX = Math.clamp(x + xx, 0, img.width - 1);
                        final int clampedY = Math.clamp(y + yy, 0, img.height - 1);

                        gradientsX += imgEdges.pixels[clampedX][clampedY] * kernel[yy + winMid][xx + winMid];
                    }
                }


                gradientsY = 0;

                for (int xx = -winMid; xx < winMid + 1; xx++) {
                    for (int yy = -winMid; yy < winMid + 1; yy++) {
                        final int clampedX = Math.clamp(x + xx, 0, img.width- 1);
                        final int clampedY = Math.clamp(y + yy, 0, img.height - 1);

                        gradientsY += imgEdges.pixels[clampedX][clampedY] * kernel[xx + winMid][yy + winMid];
                    }
                }

                int angleDeg = (int) Math.toDegrees(Math.atan2(gradientsY, gradientsX));

                if (angleDeg < 0) {
                    angleDeg += 180;
                }

                gradientsAngleDeg[x][y] = angleDeg;
                magnitude[x][y] = (int) Math.sqrt(gradientsX * gradientsX + gradientsY * gradientsY);
            }
        }

        imgEdges = new Image(img.depth, img.width, img.height);

        for (int x = 1; x < img.width - 1; x++) {
            for (int y = 1; y < img.height - 1; y++) {
                int neighbour1;
                int neighbour2;

                if (gradientsAngleDeg[x][y] < 22.5 || gradientsAngleDeg[x][y] > 157.5) {
                    neighbour1 = magnitude[x][y - 1];
                    neighbour2 = magnitude[x][y + 1];
                } else if (gradientsAngleDeg[x][y] > 112.5) {
                    neighbour1 = magnitude[x - 1][y + 1];
                    neighbour2 = magnitude[x + 1][y - 1];
                } else if (gradientsAngleDeg[x][y] < 67.5) {
                    neighbour1 = magnitude[x - 1][y - 1];
                    neighbour2 = magnitude[x + 1][y + 1];
                } else {
                    neighbour1 = magnitude[x - 1][y];
                    neighbour2 = magnitude[x + 1][y];
                }

                if (magnitude[x][y] >= neighbour1 && magnitude[x][y] >= neighbour2) {
                    imgEdges.pixels[x][y] = magnitude[x][y];
                }
            }
        }

        rescaleImg(imgEdges);

        for (int x = 1; x < img.width - 1; x++) {
            for (int y = 1; y < img.height - 1; y++) {
                if (imgEdges.pixels[x][y] < imgEdges.depth / 8) {
                    imgEdges.pixels[x][y] = 0;
                }
            }
        }

        return imgEdges;
    }

    private static Image sobelEdgeDetector(Image img) {
        final Image imgEdges = new Image(img.depth, img.width, img.height);

        final int winSize = 3;
        final int winMid = winSize / 2;
        final int[][] kernel = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};

        int horizontal;
        int vertical;

        for (int x = 0; x < img.width; x++) {
            for (int y = 0; y < img.height; y++) {
                horizontal = 0;

                for (int xx = -winMid; xx < winMid + 1; xx++) {
                    for (int yy = -winMid; yy < winMid + 1; yy++) {
                        final int clampedX = Math.clamp(x + xx, 0, img.width - 1);
                        final int clampedY = Math.clamp(y + yy, 0, img.height - 1);

                        horizontal += img.pixels[clampedX][clampedY] * kernel[yy + winMid][xx + winMid];
                    }
                }

                vertical = 0;

                for (int xx = -winMid; xx < winMid + 1; xx++) {
                    for (int yy = -winMid; yy < winMid + 1; yy++) {
                        final int clampedX = Math.clamp(x + xx, 0, img.width- 1);
                        final int clampedY = Math.clamp(y + yy, 0, img.height - 1);

                        vertical += img.pixels[clampedX][clampedY] * kernel[xx + winMid][yy + winMid];
                    }
                }

                imgEdges.pixels[x][y] = (int) Math.sqrt(horizontal * horizontal + vertical * vertical);
            }
        }

        return imgEdges;
    }

    private static void rescaleImg(Image img) {
        int min = img.depth;
        int max = 0;

        for (int x = 0; x < img.width; x++) {
            for (int y = 0; y < img.height; y++) {
                min = Math.min(img.pixels[x][y], min);
                max = Math.max(img.pixels[x][y], max);
            }
        }

        if (max - min == 0) {
            System.out.println("WARN: attempted to rescale image but range of values was 0");
            return;
        }

        for (int x = 0; x < img.width; x++) {
            for (int y = 0; y < img.height; y++) {
                img.pixels[x][y] = (img.pixels[x][y] - min) * img.depth / (max - min);
            }
        }
    }

    private static Image gaussianBlur(Image imgToBlur, double sigma) {
        final Image imgBlurred = new Image(imgToBlur.depth, imgToBlur.width, imgToBlur.height);
        final double[][] imgTemp = new double[imgToBlur.width][imgToBlur.height];
        final double[] kernel = generate1dGaussianKernel(sigma);

        for (int x = 0; x < imgToBlur.width; x++) {
            for (int y = 0; y < imgToBlur.height; y++) {
                for (int xx = -gaussianWinMid; xx < gaussianWinMid + 1; xx++) {
                    final int clampedX = Math.clamp(x + xx, 0, imgToBlur.width - 1);

                    imgTemp[x][y] += imgToBlur.pixels[clampedX][y] * kernel[xx + gaussianWinMid];
                }
            }
        }

        double sum;

        for (int x = 0; x < imgToBlur.width; x++) {
            for (int y = 0; y < imgToBlur.height; y++) {
                sum = 0;

                for (int yy = -gaussianWinMid; yy < gaussianWinMid + 1; yy++) {
                    final int clampedY = Math.clamp(y + yy, 0, imgToBlur.height - 1);

                    sum += imgTemp[x][clampedY] * kernel[yy + gaussianWinMid];
                }

                imgBlurred.pixels[x][y] = (int) sum;
            }
        }

        System.out.println("INFO: applied gaussian blur");

        return imgBlurred;
    }

    private static double[] generate1dGaussianKernel(double sigma) {
        final double denominator = 2.0f * sigma * sigma;
        final double[] kernel = new double[gaussianWinSize];

        double sum = 0;

        for (int x = -gaussianWinMid; x < gaussianWinMid + 1; x++) {
            kernel[x + gaussianWinMid] = Math.exp(-(x * x) / denominator);
            sum += kernel[x + gaussianWinMid];
        }

        for (int x = 0; x < gaussianWinSize; x++) {
            kernel[x] /= sum;
        }

        System.out.println("INFO: generated gaussian kernel with sigma " + sigma);
        System.out.println(Arrays.toString(kernel));

        return kernel;
    }

    private static double[][] houghTransformLine(Image img) {
        final int imgWidthHalf = img.width / 2;
        final int imgHeightHalf = img.height / 2;
        final int lineAccumulatorWidth = (int) Math.ceil(Math.sqrt(img.width * img.width + img.height * img.height));
        final int lineAccumulatorWidthHalf = lineAccumulatorWidth / 2;
        final int lineAccumulatorHeight = 180;

        final double[][] lineAccumulator = new double[lineAccumulatorWidth][lineAccumulatorHeight];

        for (int x = momentsWinMid; x < img.width - momentsWinMid; x++) {
            for (int y = momentsWinMid; y < img.height - momentsWinMid; y++) {
                if (img.pixels[x][y] == 0) {
                    continue;
                }

                // Math.PI / 2 is offset for estimating theta for polar coordinates, since its perpendicular to
                // moments orientation
                double estimatedThetaRad = getEstimatedThetaRad(img, x, y) + Math.PI / 2;

                int lower = (int) (Math.toDegrees(estimatedThetaRad) - deltaThetaDeg);
                int upper = (int) (Math.toDegrees(estimatedThetaRad) + deltaThetaDeg) + 1;

                final int relativeX = x - imgWidthHalf;
                final int relativeY = y - imgHeightHalf;

                for (int thetaDeg = lower; thetaDeg <= upper; thetaDeg++) {
                    int rho;

                    if (thetaDeg < 0) {
                        rho = (int) (relativeX * cosDeg[thetaDeg + 360] + relativeY * sinDeg[thetaDeg + 360]);
                    } else {
                        rho = (int) (relativeX * cosDeg[thetaDeg] + relativeY * sinDeg[thetaDeg]);
                    }

                    int adjustedThetaDeg = thetaDeg;

                    if (adjustedThetaDeg < 0) {
                        adjustedThetaDeg += 180;
                        rho *= -1;
                    } else if (adjustedThetaDeg >= 180) {
                        adjustedThetaDeg -= 180;
                        rho *= -1;
                    }

                    lineAccumulator[rho + lineAccumulatorWidthHalf][adjustedThetaDeg] += img.pixels[x][y];
                }
            }
        }

        System.out.println("INFO: generated line accumulator space with size: " + lineAccumulator.length + "x" + lineAccumulator[0].length);

        return lineAccumulator;
    }

    private static double getEstimatedThetaRad(Image img, int x, int y) {
        double m00 = 0.0f;
        double m10 = 0.0f;
        double m01 = 0.0f;

        for (int xx = -momentsWinMid; xx < momentsWinMid + 1; xx++) {
            for (int yy = -momentsWinMid; yy < momentsWinMid + 1; yy++) {
                final int weight = img.pixels[x + xx][y + yy];

                m00 += weight;
                m10 += xx * weight;
                m01 += yy * weight;
            }
        }

        final double centroidX = m10 / m00;
        final double centroidY = m01 / m00;

        double mu11 = 0.0f;
        double mu20 = 0.0f;
        double mu02 = 0.0f;

        for (int xx = -momentsWinMid; xx < momentsWinMid + 1; xx++) {
            for (int yy = -momentsWinMid; yy < momentsWinMid + 1; yy++) {
                final int weight = img.pixels[x + xx][y + yy];

                mu11 += (xx - centroidX) * (yy - centroidY) * weight;
                mu20 += (xx - centroidX) * (xx - centroidX) * weight;
                mu02 += (yy - centroidY) * (yy - centroidY) * weight;
            }
        }

        return 0.5 * Math.atan2(2 * mu11, mu20 - mu02);
    }

    private static Image getImgLineAccumulator(double[][] lineAccumulator) {
        final int lineAccumulatorImgDepth = 255;

        Image image = new Image(lineAccumulatorImgDepth, lineAccumulator.length, lineAccumulator[0].length);

        double min = 0;
        double max = 0;

        for (double[] radii : lineAccumulator) {
            for (double value : radii) {
                min = Math.min(min, value);
                max = Math.max(max, value);
            }
        }

        if (min == max) {
            System.out.println("WARN: attempted to rescale image but range of values was 0");
            return image;
        }

        for (int x = 0; x < image.width; x++) {
            for (int y = 0; y < image.height; y++) {
                image.pixels[x][y] = (int) ((lineAccumulator[x][y] - min) * image.depth / (max - min));
            }
        }

        return image;
    }

    private static void getPeaksLineAccumulator(double[][] lineAccumulator, ArrayList<Line> lines) {
        final int lineAccumulatorWidthHalf = lineAccumulator.length / 2;

        double maxVotes = 0;

        for (double[] radii : lineAccumulator) {
            for (double value : radii) {
                maxVotes = Math.max(maxVotes, value);
            }
        }

        final int threshold = (int) (linePeakThreshold * maxVotes);

        System.out.println("INFO: line accumulator max votes - " + maxVotes + ", threshold - " + threshold);

        boolean isPeak;
        double value;

        for (int rho = 0; rho < lineAccumulator.length; rho++) {
            for (int thetaDeg = 0; thetaDeg < lineAccumulator[0].length; thetaDeg++) {
                value = lineAccumulator[rho][thetaDeg];

                isPeak = true;

                checkIsPeak:
                for (int xx = -peaksWinMid; xx < peaksWinMid + 1; xx++) {
                    for (int yy = -peaksWinMid; yy < peaksWinMid + 1; yy++) {
                        final int checkX = Math.clamp(rho + xx, 0, lineAccumulator.length - 1);
                        final int checkY = Math.clamp(thetaDeg + yy, 0, lineAccumulator[0].length - 1);

                        if (lineAccumulator[checkX][checkY] > value) {
                            isPeak = false;
                            break checkIsPeak;
                        }
                    }
                }

                if (isPeak && value > threshold) {
                    lines.add(new Line(rho - lineAccumulatorWidthHalf, thetaDeg));
                }
            }
        }

        System.out.println("INFO: number of line peaks found: " + lines.size());
    }

    private static ImagePPM overlayLinesFromPeaks(Image img, ArrayList<Line> lines) {
        final double imgWidthHalf = (double) img.width / 2;
        final double imgHeightHalf = (double) img.height / 2;

        ImagePPM imgOverlay = new ImagePPM(img.depth, img.width, img.height);

        for (int c = 0; c < 3; c++) {
            for (int x = 0; x < img.width; x++) {
                System.arraycopy(img.pixels[x], 0, imgOverlay.pixels[c][x], 0, img.height);
            }
        }

        for (Line line : lines) {
            final double rho = line.rho;
            final double thetaRad = Math.toRadians(line.thetaDeg);

            double gradient = -1 / Math.tan(thetaRad);

            if (Math.abs(Math.tan(thetaRad)) < 0.0001) {
                gradient = 10000;
            }

            final double xCoordinate = imgWidthHalf + rho * cosDeg[line.thetaDeg];
            final double yCoordinate = imgHeightHalf + rho * sinDeg[line.thetaDeg];
            final double yIntercept = yCoordinate - gradient * xCoordinate;

            if (Math.abs(gradient) > 1) {
                for (int y = 0; y < imgOverlay.height; y++) {
                    final int x = (int) ((y - yIntercept) / gradient);

                    if (x < 0 || x >= imgOverlay.width) {
                        continue;
                    }

                    imgOverlay.pixels[0][x][y] = 0;
                    imgOverlay.pixels[1][x][y] = imgOverlay.depth;
                    imgOverlay.pixels[2][x][y] = 0;
                }
            } else {
                for (int x = 0; x < imgOverlay.width; x++) {
                    final int y = (int) (gradient * x + yIntercept);

                    if (y < 0 || y >= imgOverlay.height) {
                        continue;
                    }

                    imgOverlay.pixels[0][x][y] = 0;
                    imgOverlay.pixels[1][x][y] = imgOverlay.depth;
                    imgOverlay.pixels[2][x][y] = 0;
                }
            }
        }

        return imgOverlay;
    }

    private static double[][][] houghTransformSquare(double[][] lineAccumulator) {
        final int lineAccumulatorHalfWidth = lineAccumulator.length / 2;
        final int squareAccumulatorWidth = (int) Math.sqrt((double) lineAccumulator.length * lineAccumulator.length / 2);
        final int squareAccumulatorHalfWidth = squareAccumulatorWidth / 2;
        final int squareAccumulatorDepth = 90;

        final double[][][] squareAccumulator = new double[squareAccumulatorWidth][squareAccumulatorWidth][squareAccumulatorDepth];

        for (int x = 0; x < squareAccumulator.length; x++) {
            for (int y = 0; y < squareAccumulator[0].length; y++) {
                for (int thetaDeg = 0; thetaDeg < squareAccumulator[0][0].length; thetaDeg++) {
                    final int relativeSquareCentreX = x - squareAccumulatorHalfWidth;
                    final int relativeSquareCentreY = y - squareAccumulatorHalfWidth;

                    final double[] squareSidesX = {squareSizeHalf, 0.0f, -squareSizeHalf, 0.0f};
                    final double[] squareSidesY = {0.0f, squareSizeHalf, 0.0f, -squareSizeHalf};
                    final int[] squareSidesThetaDeg = {thetaDeg, thetaDeg + 90, thetaDeg + 180, thetaDeg + 270};
                    final int[] squareSidesRho = {0, 0, 0, 0};
                    final double[] mins = {0, 0, 0, 0};

                    for (int i = 0; i < squareSidesX.length; i++) {
                        final double temp = squareSidesX[i] * cosDeg[thetaDeg] - squareSidesY[i] * sinDeg[thetaDeg] + relativeSquareCentreX;
                        squareSidesY[i] = squareSidesX[i] * sinDeg[thetaDeg] + squareSidesY[i] * cosDeg[thetaDeg] + relativeSquareCentreY;
                        squareSidesX[i] = temp;

                        squareSidesRho[i] = (int) (Math.abs(squareSidesX[i]) * cosDeg[thetaDeg] + Math.abs(squareSidesY[i]) * sinDeg[thetaDeg]);

                        int adjustedThetaDeg = squareSidesThetaDeg[i];

                        if (adjustedThetaDeg < 0) {
                            adjustedThetaDeg += 180;
                            squareSidesRho[i] *= -1;
                        } else if (adjustedThetaDeg >= 180) {
                            adjustedThetaDeg -= 180;
                            squareSidesRho[i] *= -1;
                        }

                        final int rho = squareSidesRho[i] + lineAccumulatorHalfWidth;

                        if (rho < 0 || rho >= lineAccumulator.length) {
                            continue;
                        }

                        mins[i] = lineAccumulator[rho][adjustedThetaDeg];
                    }

                    Arrays.sort(mins);

                    squareAccumulator[x][y][thetaDeg] = mins[0];
                }
            }
        }

        System.out.println("INFO: generated square accumulator space with size: " + squareAccumulator.length + "x" + squareAccumulator[0].length + "x" + squareAccumulator[0][0].length);

        return squareAccumulator;
    }

    private static void getPeaksSquareAccumulator(double[][][] squareAccumulator, ArrayList<Square> squares) {
        final int squareAccumulatorWidthHalf = squareAccumulator.length / 2;

        double maxVotes = 0;

        for (double[][] intsX : squareAccumulator) {
            for (double[] intsY : intsX) {
                for (double votes : intsY) {
                    maxVotes = Math.max(maxVotes, votes);
                }
            }
        }

        final int threshold = (int) (squarePeakThreshold * maxVotes);

        System.out.println("INFO: square accumulator max votes - " + maxVotes + ", threshold - " + threshold);

        boolean isPeak;
        double value;

        for (int x = 0; x < squareAccumulator.length; x++) {
            for (int y = 0; y < squareAccumulator[0].length; y++) {
                for (int thetaDegrees = 0; thetaDegrees < squareAccumulator[0][0].length; thetaDegrees++) {
                    value = squareAccumulator[x][y][thetaDegrees];

                    if (value == 0.0f) {
                        continue;
                    }

                    isPeak = true;

                    checkIsPeak:
                    for (int xx = -peaksWinMid; xx < peaksWinMid + 1; xx++) {
                        for (int yy = -peaksWinMid; yy < peaksWinMid + 1; yy++) {
                            for (int tt = -peaksWinMid; tt < peaksWinMid + 1; tt++) {
                                final int checkX = x + xx;
                                final int checkY = y + yy;
                                final int checkThetaDegrees = thetaDegrees + tt;

                                if (checkX < 0 || checkX >= squareAccumulator.length || checkY < 0 || checkY >= squareAccumulator[0].length || checkThetaDegrees < 0 || checkThetaDegrees >= squareAccumulator[0][0].length) {
                                    continue;
                                }

                                if (squareAccumulator[checkX][checkY][checkThetaDegrees] > value) {
                                    isPeak = false;
                                    break checkIsPeak;
                                }
                            }
                        }
                    }

                    if (isPeak && value > threshold) {
                        squares.add(new Square(x - squareAccumulatorWidthHalf, y - squareAccumulatorWidthHalf, thetaDegrees));
                    }
                }
            }
        }

        System.out.println("INFO: number of square peaks found: " + squares.size());
    }

    private static ArrayList<SquarePixels> getSquaresPixelsFromPeaks(Image image, ArrayList<Square> squares) {
        final ArrayList<SquarePixels> squarePixelsList = new ArrayList<>();

        for (Square square : squares) {
            final int relativeSquareCentreX = square.centreX + image.width / 2;
            final int relativeSquareCentreY = square.centreY + image.height / 2;
            final double squareRotationRadians = Math.toRadians(square.rotationDeg);

            if (relativeSquareCentreX < 0 || relativeSquareCentreX >= image.width || relativeSquareCentreY < 0 || relativeSquareCentreY >= image.height) {
                continue;
            }

            final double[] squareSidesX = {squareSizeHalf, -squareSizeHalf, -squareSizeHalf, squareSizeHalf};
            final double[] squareSidesY = {squareSizeHalf, squareSizeHalf, -squareSizeHalf, -squareSizeHalf};
            final int[] squareSidesXInt = {0, 0, 0, 0};
            final int[] squareSidesYInt = {0, 0, 0, 0};

            for (int j = 0; j < squareSidesX.length; j++) {
                final double temp = squareSidesX[j] * Math.cos(squareRotationRadians) - squareSidesY[j] * Math.sin(squareRotationRadians);
                squareSidesY[j] = squareSidesX[j] * Math.sin(squareRotationRadians) + squareSidesY[j] * Math.cos(squareRotationRadians);
                squareSidesX[j] = temp;

                squareSidesXInt[j] = (int) squareSidesX[j] + relativeSquareCentreX;
                squareSidesYInt[j] = (int) squareSidesY[j] + relativeSquareCentreY;
            }

            SquarePixels squarePixels = getSquarePixels(image, squareSidesXInt, squareSidesYInt);

            squarePixelsList.add(squarePixels);
        }

        return squarePixelsList;
    }

    private static SquarePixels getSquarePixels(Image image, int[] squareSidesXInt, int[] squareSidesYInt) {
        SquarePixels squarePixels = new SquarePixels();

        for (int j = 0; j < squareSidesXInt.length; j++) {
            int startX = squareSidesXInt[j];
            int startY = squareSidesYInt[j];

            int checkIndex = j + 1;

            if (checkIndex >= squareSidesXInt.length) {
                checkIndex -= squareSidesXInt.length;
            }

            int endX = squareSidesXInt[checkIndex];
            int endY = squareSidesYInt[checkIndex];

            double gradient = (double) (endY - startY) / (endX - startX);
            final double yIntercept = endY - gradient * endX;

            if (endX == startX) {
                gradient = 10000;
            }

            if (Math.abs(gradient) > 1) {
                final int y0 = Math.min(startY, endY);
                final int y1 = Math.max(startY, endY);

                for (int y = y0; y <= y1; y++) {
                    final int x = (int) ((y - yIntercept) / gradient);

                    if (x < 0 || x >= image.width || y < 0 || y >= image.height) {
                        continue;
                    }

                    squarePixels.add(x, y);
                }
            } else {
                final int x0 = Math.min(startX, endX);
                final int x1 = Math.max(startX, endX);

                for (int x = x0; x <= x1; x++) {
                    final int y = (int) (gradient * x + yIntercept);

                    if (x < 0 || x >= image.width || y < 0 || y >= image.height) {
                        continue;
                    }

                    squarePixels.add(x, y);
                }
            }
        }
        return squarePixels;
    }

    private static ImagePPM overlaySquaresFromPixels(Image img, ArrayList<SquarePixels> squarePixelsList) {
        final ImagePPM imgOverlay = new ImagePPM(img.depth, img.width, img.height);

        for (int c = 0; c < 3; c++) {
            for (int x = 0; x < img.width; x++) {
                System.arraycopy(img.pixels[x], 0, imgOverlay.pixels[c][x], 0, img.height);
            }
        }

        for (SquarePixels squarePixels : squarePixelsList) {
            for (int i = 0; i < squarePixels.x.size(); i++) {
                final int x = squarePixels.x.get(i);
                final int y = squarePixels.y.get(i);

                imgOverlay.pixels[0][x][y] = img.depth;
                imgOverlay.pixels[1][x][y] = 0;
                imgOverlay.pixels[2][x][y] = 0;
            }
        }

        return imgOverlay;
    }

    private static ImagePPM overlaySquaresFromPixelsWithBackProj(Image img, Image dog, ArrayList<SquarePixels> squarePixelsList) {
        final ImagePPM imgOverlay = new ImagePPM(img.depth, img.width, img.height);

        for (int c = 0; c < 3; c++) {
            for (int x = 0; x < img.width; x++) {
                System.arraycopy(img.pixels[x], 0, imgOverlay.pixels[c][x], 0, img.height);
            }
        }

        int maxVotes = 0;
        int sum;

        for (SquarePixels squarePixels : squarePixelsList) {
            sum = 0;

            for (int i = 0; i < squarePixels.x.size(); i++) {
                sum += dog.pixels[squarePixels.x.get(i)][squarePixels.y.get(i)];
            }

            maxVotes = Math.max(maxVotes, sum);
        }

        final int threshold = (int) (maxVotes * backProjThreshold);

        System.out.println("INFO: square accumulator back projection max votes - " + maxVotes + ", threshold - " + threshold);

        int count = 0;

        for (SquarePixels squarePixels : squarePixelsList) {
            sum = 0;

            for (int j = 0; j < squarePixels.x.size(); j++) {
                sum += dog.pixels[squarePixels.x.get(j)][squarePixels.y.get(j)];
            }

            if (sum < threshold) {
                continue;
            }

            count++;

            for (int j = 0; j < squarePixels.x.size(); j++) {
                final int x = squarePixels.x.get(j);
                final int y = squarePixels.y.get(j);

                imgOverlay.pixels[0][x][y] = 0;
                imgOverlay.pixels[1][x][y] = 0;
                imgOverlay.pixels[2][x][y] = img.depth;
            }
        }

        System.out.println("INFO: squares passing back projection found: " + count);

        return imgOverlay;
    }
}

class Line {
    public int rho;
    public int thetaDeg;

    Line (int rho, int thetaDeg) {
        this.rho = rho;
        this.thetaDeg = thetaDeg;
    }
}

class Square {
    public int centreX;
    public int centreY;
    public int rotationDeg;

    Square(int centreX, int centreY, int rotationDeg) {
        this.centreX = centreX;
        this.centreY = centreY;
        this.rotationDeg = rotationDeg;
    }
}

class SquarePixels {
    public ArrayList<Integer> x = new ArrayList<>();
    public ArrayList<Integer> y = new ArrayList<>();

    public void add(int x, int y) {
        this.x.add(x);
        this.y.add(y);
    }
}