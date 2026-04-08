# Square Detection
Implementation to detect squares in PGM images using hough transforms written in java.

- Edge detection is performed on the image (difference of guassian, sobel)
- For each edge pixel a moment is calculated to get a direction
- This direction and intensity of pixel is used in a hough transform for lines
- Line accumulator is drawn for a visual guide
- Another hough transform is performed using the lines previously calculated to find squares
- Square accumulator is drawn for a visual guide
