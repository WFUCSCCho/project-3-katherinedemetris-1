/******************************************************************************
 * @file: HousingPricesData.java
 * @description: This class represents a single housing price data entry with
 *               various attributes such as suburb, address, price, etc. It
 *               implements Comparable for sorting based on suburb names.
 * @author: Katherine Demetris
 * @date: November 10, 2024
 ******************************************************************************/

public class HousingPricesData implements Comparable<HousingPricesData> {
    private String suburb;
    private String address;
    private Integer rooms;
    private Character type;
    private Integer price;
    private Character method;
    private String sellerG;
    private String date;
    private Integer postcode;
    private String regionName;
    private Integer  propertyCount;
    private Double  distance;
    private String councilArea;

    // Default constructor
    public HousingPricesData() {
        suburb = " ";
        address = " ";
        rooms = 0;
        type = ' ';
        price = 0;
        method = ' ';
        sellerG = " ";
        date = " ";
        postcode = 0;
        regionName = " ";
        propertyCount = 0;
        distance = 0.0;
        councilArea = " ";
    }

    // Parametrized constructor
    public HousingPricesData(String suburb, String address, int rooms, char type, int price, char method, String sellerG, String date, int postcode, String regionName, int propertyCount, double distance, String councilArea) {
        this.suburb = suburb != null ? suburb : "";
        this.address = address != null ? address : "";
        this.rooms = rooms;
        this.type = type;
        this.price = price;
        this.method = method;
        this.sellerG = sellerG != null ? sellerG : "";
        this.date = date != null ? date : "";
        this.postcode = postcode;
        this.regionName = regionName != null ? regionName : "";
        this.propertyCount = propertyCount;
        this.distance = distance;
        this.councilArea = councilArea != null ? councilArea : "";
    }

    // Copy constructor
    public HousingPricesData(HousingPricesData other) {
        this.suburb = other.suburb;
        this.address = other.address;
        this.rooms = other.rooms;
        this.type = other.type;
        this.price = other.price;
        this.method = other.method;
        this.sellerG = other.sellerG;
        this.date = other.date;
        this.postcode = other.postcode;
        this.regionName = other.regionName;
        this.propertyCount = other.propertyCount;
        this.distance = other.distance;
        this.councilArea = other.councilArea;
    }


    //TODO: toString method --> this method should return a string representation of the object in the format of your choosing.
    @Override
    public String toString() {
        return String.format(
                "%s: Address: %s, Rooms: %d, Type: %c, Price: $%,d, Method: %c, " +
                        "Seller: %s, Date: %s, Postcode: %d, Region: %s, " +
                        "Property Count: %d, Distance: %.1f km, Council Area: %s\n",
                suburb, address, rooms, type, price, method, sellerG, date,
                postcode, regionName, propertyCount, distance, councilArea
        );
    }

    // This method should return true if the object is equal to the object passed as a parameter.
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null || getClass() != obj.getClass()) {
            return false;
        }
        HousingPricesData house = (HousingPricesData) obj;
        return rooms.equals(house.rooms) && price.equals(house.price) &&
                suburb.equals(house.suburb) && address.equals(house.address) &&
                type.equals(house.type) && method.equals(house.method) &&
                sellerG.equals(house.sellerG) && date.equals(house.date) &&
                postcode.equals(house.postcode) && regionName.equals(house.regionName) &&
                propertyCount.equals(house.propertyCount) && distance.equals(house.distance) &&
                councilArea.equals(house.councilArea);
    }


    //TODO: Comparable interface for HousingPrices data
    @Override
    public int compareTo(HousingPricesData other) {
        // Compare by price from (highest to lowest)
        return -1 * Integer.compare(this.price, other.price);
    }

    // Getters and Setters
    public String getSuburb() {return suburb;}
    public void setSuburb(String suburb) {this.suburb = suburb;}

    public String getAddress() {return address;}
    public void setAddress(String address) {this.address = address;}

    public int getRooms() {return rooms;}
    public void setRooms(int rooms) {this.rooms = rooms;}

    public char getType() {return type;}
    public void setType(char type) {this.type = type;}

    public int getPrice() {return price;}
    public void setPrice(int price) {this.price = price;}

    public char getMethod() {return method;}
    public void setMethod(char method) {this.method = method;}

    public String getSellerG() {return sellerG;}
    public void setSellerG(String sellerG) {this.sellerG = sellerG;}

    public String getDate() {return date;}
    public void setDate(String date) {this.date = date;}

    public int getPostcode() {return postcode;}
    public void setPostcode(int postcode) {this.postcode = postcode;}

    public String getRegionName() {return regionName;}
    public void setRegionName(String regionName) {this.regionName = regionName;}

    public int getPropertyCount() {return propertyCount;}
    public void setPropertyCount(int propertyCount) {this.propertyCount = propertyCount;}

    public double getDistance() {return distance;}
    public void setDistance(double distance) {this.distance = distance;}

    public String getCouncilArea() {return councilArea;}
    public void setCouncilArea(String councilArea) {this.councilArea = councilArea;}


}


