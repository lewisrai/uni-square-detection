/* rather simplistic implementation of PPM reader/writer
 *
 * Paul Rosin, 2003
 *
 * update to use 2 input streams - avoiding deprecation errors
 * Matt Gee, 2004
 *
 * constructor added to allow image size to be specified
 * Graham Daniell, 2004
 *
 * pixels array overwritten during read to match image dimensions
 * Jin Huang, 2022
 *
 * check magic number
 * Simeon Alcolado, 2023
 *
 * fix bug with magic number code & added support for PPM files without comments
 * Constantinos Psomadakis, 2024
 */

 import java.io.*;
 import java.util.*;

 public class ImagePPM
 {
     public int [][][] pixels;
     public int depth,width,height;

     public ImagePPM()
     {
         pixels = new int[3][1500][1500];
         depth = width = height = 0;
     }

     public ImagePPM(int inDepth, int inWidth, int inHeight)
     {
         pixels = new int[3][inWidth][inHeight];
         width = inWidth;
         height = inHeight;
         depth = inDepth;
     }

     public void ReadPPM(String fileName)
     {
         String line;
         StringTokenizer st;

         try {
             BufferedReader in =
               new BufferedReader(new InputStreamReader(
                 new BufferedInputStream(
                   new FileInputStream(fileName))));

             DataInputStream in2 =
               new DataInputStream(
                 new BufferedInputStream(
                   new FileInputStream(fileName)));

             // read PPM image header
             line = in.readLine();
             in2.skip((line+"\n").getBytes().length);
             if (!Objects.equals(line,"P6")) {
                 System.out.println("ERROR: image is not correct PPM format");
                 System.exit(0);
             }

             // skip comments
             do {
                 line = in.readLine();
                 in2.skip((line+"\n").getBytes().length);
             } while (line.charAt(0) == '#');

             // the current line has dimensions
             st = new StringTokenizer(line);
             width = Integer.parseInt(st.nextToken());
             height = Integer.parseInt(st.nextToken());

             // overwrite pixels array to the correct dimensions
             pixels = new int[3][width][height];

             // next line has pixel depth
             line = in.readLine();
             in2.skip((line+"\n").getBytes().length);
             st = new StringTokenizer(line);
             depth = Integer.parseInt(st.nextToken());
             if (depth != 255) {
                 System.out.println("Error: depth = "+depth+" (instead of 255)");
                 System.exit(0);
             }

             // read pixels now
             for (int y = 0; y < height; y++)
                 for (int x = 0; x < width; x++)
                     for (int i = 0; i < 3; i++)
                         pixels[i][x][y] = in2.readUnsignedByte();

             in.close();
             in2.close();
         } catch(ArrayIndexOutOfBoundsException e) {
             System.out.println("Error: image in "+fileName+" too big");
         } catch(FileNotFoundException e) {
             System.out.println("Error: file "+fileName+" not found");
         } catch(IOException e) {
             System.out.println("Error: end of stream encountered when reading "+fileName);
         }
     }

     public void WritePPM(String fileName)
     {
         String line;
         StringTokenizer st;

         try {
             DataOutputStream out =
               new DataOutputStream(
                 new BufferedOutputStream(
                   new FileOutputStream(fileName)));

             out.writeBytes("P6\n");
             out.writeBytes("#created by Paul Rosin\n");
             out.writeBytes(width+" "+height+"\n255\n");

             for (int y = 0; y < height; y++)
                 for (int x = 0; x < width; x++)
                     for (int i = 0; i < 3; i++)
                         out.writeByte((byte)pixels[i][x][y]);

             out.close();
         } catch(IOException e) {
             System.out.println("ERROR: cannot write output file");
         }
     }
 }

