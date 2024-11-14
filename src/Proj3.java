/******************************************************************************
 * @file: Proj3.java
 * @description: This program implements various sorting algorithms to analyze
 *              performance on housing price data
 * @author: Katherine Demetris
 * @date: November 14, 2024
 ******************************************************************************/

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class Proj3 {

    // Add these static variables at the class level
    private static int mergeSortComparisons = 0;
    private static int quickSortComparisons = 0;
    private static int heapSortComparisons = 0;

    // Sorting Method declarations

    /**
     * MERGESORT
     * Internal method that makes recursive calls.
     *
     * @param a     an ArrayList of Comparable items.
     * @param left  the left-most index of the subarray.
     * @param right the right-most index of the subarray.
     */
    public static <T extends Comparable> void mergeSort(ArrayList<T> a, int left, int right) {
        // Finish Me
        if (left < right) {
            // Find the midpoint of the array
            int mid = (left + right) / 2;

            // Recursively sort the left half
            mergeSort(a, left, mid);

            // Recursively sort the right half
            mergeSort(a, mid + 1, right);

            // Merge the sorted halves
            merge(a, left, mid, right);
        }
    }

    /**
     * MERGE
     * Internal method that merges two sorted halves of a subarray.
     *
     * @param a     an ArrayList of Comparable items.
     * @param left  the left-most index of the subarray.
     * @param mid   the index of the start of the second half.
     * @param right the right-most index of the subarray.
     */
    public static <T extends Comparable> void merge(ArrayList<T> a, int left, int mid, int right) {
        // Finish Me
        ArrayList<T> leftArray = new ArrayList<>();
        ArrayList<T> rightArray = new ArrayList<>();

        for (int i = left; i <= mid; i++) {
            leftArray.add(a.get(i));
        }
        for (int j = mid + 1; j <= right; j++) {
            rightArray.add(a.get(j));
        }

        int i = 0, j = 0, k = left;

        while (i < leftArray.size() && j < rightArray.size()) {
            mergeSortComparisons++; // Add counter here
            if (leftArray.get(i).compareTo(rightArray.get(j)) <= 0) {
                a.set(k++, leftArray.get(i++));
            } else {
                a.set(k++, rightArray.get(j++));
            }
        }

        // Copy remaining elements
        while (i < leftArray.size()) {
            a.set(k++, leftArray.get(i++));
        }
        while (j < rightArray.size()) {
            a.set(k++, rightArray.get(j++));
        }
    }

    /**
     * QUICKSORT
     * Internal quicksort method that makes recursive calls.
     * Uses median-of-three partitioning and a cutoff of 10.
     *
     * @param a     an arrayList of Comparable items.
     * @param left  the left-most index of the subarray.
     * @param right the right-most index of the subarray.
     */
    public static <T extends Comparable> void quickSort(ArrayList<T> a, int left, int right) {
        // Finish Me
        if (left + 10 <= right) {
            int mid = (left + right) / 2;

            quickSortComparisons++;
            if (a.get(mid).compareTo(a.get(left)) < 0)
                swap(a, left, mid);
            quickSortComparisons++;
            if (a.get(right).compareTo(a.get(left)) < 0)
                swap(a, left, right);
            quickSortComparisons++;
            if (a.get(right).compareTo(a.get(mid)) < 0)
                swap(a, mid, right);

            swap(a, mid, right - 1);
            T pivot = a.get(right - 1);

            int i = left;
            int j = right - 1;

            while (true) {
                while (i < right - 1) {
                    quickSortComparisons++;
                    if (a.get(++i).compareTo(pivot) >= 0) break;
                }
                while (j > left) {
                    quickSortComparisons++;
                    if (a.get(--j).compareTo(pivot) <= 0) break;
                }
                if (i >= j) break;
                swap(a, i, j);
            }

            swap(a, i, right - 1);
            quickSort(a, left, i - 1);
            quickSort(a, i + 1, right);
        } else {
            for (int i = left; i < right; i++) {
                for (int j = left; j < right - (i - left); j++) {
                    quickSortComparisons++;
                    if (a.get(j).compareTo(a.get(j + 1)) > 0) {
                        swap(a, j, j + 1);
                    }
                }
            }
        }
    }

    /**
     * PARTITION
     * Partitions the ArrayList around a pivot element.
     * Elements less than or equal to the pivot are moved to the left of the pivot's final position,
     * and elements greater than the pivot are moved to the right.
     *
     * @param a     the ArrayList to be partitioned
     * @param left  the starting index of the partitioning range
     * @param right the ending index of the partitioning range (pivot is typically chosen from here)
     * @return the index of the pivot after partitioning
     */
    public static <T extends Comparable> int partition(ArrayList<T> a, int left, int right) {
        // Finish Me
        T pivot = a.get(right); // Select the rightmost element as the pivot
        int i = left - 1;

        for (int j = left; j < right; j++) {
            if (a.get(j).compareTo(pivot) <= 0) {
                i++;
                swap(a, i, j);
            }
        }

        // Place pivot at its correct position
        swap(a, i + 1, right);
        return i + 1;
    }

    static <T> void swap(ArrayList<T> a, int i, int j) {
        T temp = a.get(i);
        a.set(i, a.get(j));
        a.set(j, temp);
    }

    /**
     * HEAPSORT
     * Standard heapsort.
     *
     * @param a an arrayList of Comparable items.
     */
    public static <T extends Comparable> void heapSort(ArrayList<T> a, int left, int right) {
        int n = a.size();
        // Build the heap
        for (int i = n / 2 - 1; i >= 0; i--) {
            heapify(a, i, n);
        }

        // Extract elements from the heap
        for (int i = n - 1; i > 0; i--) {
            swap(a, 0, i); // Move the current root to the end
            heapify(a, 0, i); // Heapify the reduced heap
        }
    }

    /**
     * HEAPIFY
     * Internal method for heapsort that is used in deleteMax and buildHeap.
     * This method ensures that the subtree rooted at index left satisfies the heap property.
     *
     * @param a     the ArrayList of Comparable items.
     * @param left  the index of the root of the subtree to heapify.
     * @param right the logical size of the heap (exclusive index).
     */
    public static <T extends Comparable> void heapify(ArrayList<T> a, int left, int right) {
        int largest = left; // Initialize largest as root
        int leftChild = 2 * left + 1; // left child
        int rightChild = 2 * left + 2; // right child

        // If left child is larger than root
        if (leftChild < right) {
            heapSortComparisons++;
            if (a.get(leftChild).compareTo(a.get(largest)) > 0)
                largest = leftChild;
        }

        // If right child is larger than largest so far
        if (rightChild < right) {
            heapSortComparisons++;
            if (a.get(rightChild).compareTo(a.get(largest)) > 0)
                largest = rightChild;
        }

        // If largest is not root
        if (largest != left) {
            swap(a, left, largest);
            heapify(a, largest, right);
        }
    }

    /**
     * BUBBLE SORT
     * Implements the Bubble Sort algorithm for sorting an ArrayList of Comparable items.
     * This algorithm compares adjacent elements and swaps them if they are in the wrong order.
     * The process repeats until no swaps are needed, indicating that the list is sorted.
     * The algorithm is optimized to stop early if the list becomes sorted before all iterations are completed.
     *
     * @param a    the ArrayList of Comparable items to be sorted.
     * @param size the number of elements in the ArrayList.
     * @return the total number of swaps made during the sorting process.
     */
    public static <T extends Comparable> int bubbleSort(ArrayList<T> a, int size) {
        //Finish me
        int comparisons = 0;
        boolean swapped;
        for (int i = 0; i < size - 1; i++) {
            swapped = false;
            for (int j = 0; j < size - i - 1; j++) {
                comparisons++; // Increment the comparison count
                if (a.get(j).compareTo(a.get(j + 1)) > 0) {
                    // Swap the elements if they're in the wrong order
                    T temp = a.get(j);
                    a.set(j, a.get(j + 1));
                    a.set(j + 1, temp);
                    swapped = true;
                }
            }
            if (!swapped) break; // If no two elements were swapped, the list is sorted
        }
        return comparisons;
    }

    /**
     * ODD-EVEN TRANSPOSITION SORT
     * Implements the Odd-Even Transposition Sort algorithm for sorting an ArrayList of Comparable items.
     * This sorting algorithm performs two phases: an "odd" phase and an "even" phase.
     * In each phase, adjacent elements are compared and swapped if they are out of order.
     * The algorithm continues until no swaps are needed, meaning the list is sorted.
     *
     * @param a the ArrayList of Comparable items to be sorted.
     * @return the total number of swaps made during the sorting process.
     */
    public static <T extends Comparable> int transpositionSort(ArrayList<T> a) {
        // Finish Me
        int size = a.size();
        int comparisons = 0;  // Change from swaps to comparisons
        boolean isSorted = false;

        // Continue sorting until the list is sorted
        while (!isSorted) {
            isSorted = true;

            // Odd phase: Compare and swap adjacent odd-indexed elements
            for (int i = 1; i < size - 1; i += 2) {
                comparisons++;  // Count the comparison
                if (a.get(i).compareTo(a.get(i + 1)) > 0) {
                    swap(a, i, i + 1);
                    isSorted = false;
                }
            }

            // Even phase: Compare and swap adjacent even-indexed elements
            for (int i = 0; i < size - 1; i += 2) {
                comparisons++;  // Count the comparison
                if (a.get(i).compareTo(a.get(i + 1)) > 0) {
                    swap(a, i, i + 1);
                    isSorted = false;
                }
            }
        }

        return comparisons;  // Return total comparisons instead of swaps
    }

    public static void main(String[] args) throws IOException {
        // Finish Me
        // Check for correct number of command-line arguments
        if (args.length != 3) {
            System.err.println("Usage: java HousingPriceSorter <input file> <sorting algorithm> <number of lines>");
            System.exit(1);
        }

        // Parse the command-line arguments
        String inputFileName = args[0];
        String sortingAlgorithm = args[1].toLowerCase();
        int numLines = Integer.parseInt(args[2]);

        // For file input
        FileInputStream inputFileNameStream = new FileInputStream(inputFileName);
        Scanner inputFileNameScanner = new Scanner(inputFileNameStream);

        // Ignore first line (header)
        inputFileNameScanner.nextLine();

        // Create ArrayList to store housing data
        ArrayList<HousingPricesData> housingPriceList = new ArrayList<>();

        // Read and parse data from the file
        int lineCount = 0;
        while (inputFileNameScanner.hasNextLine() && lineCount < numLines) {
            String line = inputFileNameScanner.nextLine();
            String[] parts = line.split(",");

            if (parts.length == 13) {
                try {
                    // Parse price by removing $ and , characters
                    String priceStr = parts[4].replaceAll("[$,]", "").trim();
                    int price = priceStr.isEmpty() ? 0 : Integer.parseInt(priceStr);

                    HousingPricesData data = new HousingPricesData(
                            parts[0], parts[1],  // suburb, address
                            parts[2].isEmpty() ? 0 : Integer.parseInt(parts[2]),  // rooms
                            parts[3].isEmpty() ? ' ' : parts[3].charAt(0),  // type
                            price,  // using cleaned price value
                            parts[5].isEmpty() ? ' ' : parts[5].charAt(0),  // method
                            parts[6], parts[7],  // sellerG, date
                            parts[8].isEmpty() ? 0 : Integer.parseInt(parts[8]),  // postcode
                            parts[9],  // regionName
                            parts[10].isEmpty() ? 0 : Integer.parseInt(parts[10]),  // propertyCount
                            parts[11].isEmpty() ? 0.0 : Double.parseDouble(parts[11]),  // distance
                            parts[12]  // councilArea
                    );
                    housingPriceList.add(data);
                    lineCount++;
                } catch (NumberFormatException e) {
                    System.err.println("Error parsing line: " + line);
                    System.err.println("Specific error: " + e.getMessage());
                }
            }
        }

        // Create sorted, randomized, and reversed copies of the data
        ArrayList<HousingPricesData> sortedData = new ArrayList<>(housingPriceList);
        ArrayList<HousingPricesData> randomizedData = new ArrayList<>(housingPriceList);
        ArrayList<HousingPricesData> reversedData = new ArrayList<>(housingPriceList);

        Collections.sort(sortedData);
        Collections.shuffle(randomizedData);
        Collections.sort(reversedData, Collections.reverseOrder());

        // Analysis file setup
        FileWriter analysisWriter;
        File analysisFile = new File("analysis.txt");
        try {
            // Check if file exists and is empty or doesn't exist
            if (!analysisFile.exists() || analysisFile.length() == 0) {
                // Create new file with header
                analysisWriter = new FileWriter("analysis.txt");
                analysisWriter.write("Sorting Method, Number of Lines, Dataset, Time (ms), Comparisons\n");
                analysisWriter.close();
            }

            // Open in append mode after ensuring header exists
            analysisWriter = new FileWriter("analysis.txt", true);

        } catch (IOException e) {
            System.err.println("Error handling analysis file: " + e.getMessage());
            throw e;
        }

        // Print header to console
        System.out.println("\nResults for sorting " + numLines + " lines:");
        System.out.println("--------------------------------------------------------");
        System.out.printf("%-8s | %-8s | %-15s | %-12s%n",
                "Algorithm", "Lines", "Time (ms)", "Comparisons");
        System.out.println("--------------------------------------------------------");

        for (String datasetType : new String[]{"Sorted", "Shuffled", "Reversed"}) {
            ArrayList<HousingPricesData> dataset = datasetType.equals("Sorted") ? sortedData :
                    datasetType.equals("Shuffled") ? randomizedData : reversedData;

            long startTime = System.nanoTime();
            int comparison = 0;

            switch (sortingAlgorithm.toLowerCase()) {
                case "bubble":
                    comparison = bubbleSort(dataset, dataset.size());
                    break;
                case "merge":
                    mergeSortComparisons = 0;  // Reset counter
                    mergeSort(dataset, 0, dataset.size() - 1);
                    comparison = mergeSortComparisons;
                    break;
                case "quick":
                    quickSortComparisons = 0;  // Reset counter
                    quickSort(dataset, 0, dataset.size() - 1);
                    comparison = quickSortComparisons;
                    break;
                case "heap":
                    heapSortComparisons = 0;  // Reset counter
                    heapSort(dataset, 0, dataset.size() - 1);
                    comparison = heapSortComparisons;
                    break;
                case "transposition":
                    comparison = transpositionSort(dataset);
                    break;
                default:
                    System.err.println("Unknown sorting algorithm: " + sortingAlgorithm);
                    System.exit(1);
            }

            long endTime = System.nanoTime();
            double elapsedTimeMs = (endTime - startTime) / 1e6;

            // Write to analysis file
            analysisWriter.append(String.format("%s, %d, %s, %.4f, %d%n",
                    sortingAlgorithm, numLines, datasetType, elapsedTimeMs, comparison));

            // After sorting, write results to sorted.txt
            FileWriter sortedWriter = new FileWriter("sorted.txt");
            sortedWriter.write("Sorted data for " + sortingAlgorithm + " algorithm:\n");
            for (HousingPricesData data : dataset) {
                sortedWriter.write(data.toString());
            }
            sortedWriter.close();

            // Print to console
            System.out.printf("%-8s | %-8d | %-15.4f | %-12d [%s]%n",
                    sortingAlgorithm, numLines, elapsedTimeMs, comparison, datasetType);
        }

        System.out.println("--------------------------------------------------------");
        System.out.println("\nResults have been written to:");
        System.out.println("- analysis.txt (appended timing data)");
        System.out.println("- sorted.txt (sorted lists)\n");

        analysisWriter.close();
    }
}

