#include <fstream>
#include <iostream>
#include <vector>
#include <string>

struct Vec2D //used for 2D positions
{
    int x, y;

    Vec2D() = default;
    Vec2D(int x, int y) : x(x), y(y) { }
    Vec2D(int dir) { //direction can be 0: LEFT, 1: UP, 2: RIGHT, 3: DOWN
            switch (dir)
            {
            case 0:
                x = -1;
                y = 0;
                break;
            case 1:
                x = 0;
                y = 1;
                break;
            case 2:
                x = 1;
                y = 0;
                break;
            case 3:
                x = 0;
                y = -1;
                break;

            default:
                x = 0;
                y = 0;
                break;
            }
    }
};

struct tuple //used to return two values in saw()
{
    int Z;
    double mean_squared;
};
inline Vec2D operator+(const Vec2D& a, const Vec2D& b) //return vector addition a + b;
{
    return {a.x + b.x, a.y + b.y};
}
inline bool operator==(const Vec2D& a, const Vec2D& b)
{
    return a.x == b.x && a.y == b.y;
}

//initialize global variables
std::ofstream out("out/saw.txt");
std::vector<Vec2D> current_path(20);
int progress_counter = 0;
unsigned long long endpoint_sum = 0;
int Z = 0;
int N_steps = 0;

bool crossing(Vec2D current_pos_candidate, int size) //returns true if current_pos_candidate crosses current_path
{
    for (int j = size-3; j > -1; j--)
    {
        if (current_pos_candidate == current_path[j])
        {
            return true;
        }
    }
    return false;
}


void saw_step(int i, int dir) //recursive step
{
    // candidate to be checked for crossing
    Vec2D current_pos_candidate = current_path[i - 1] + dir;

    // abort if candidate crosses current path
    if(crossing(current_pos_candidate, i)) return;

    // if no crossing, push candidate onto path
    current_path.emplace_back(current_pos_candidate);

    if(i == 7){ //print progress
        progress_counter++;
        std::cout << "\rProgress: " << progress_counter*100 / 2172 << " %" << std::flush;
    }

    if (i >= N_steps) //termination condition
    {

        Z++; 
        endpoint_sum += current_path[i].x*current_path[i].x + current_path[i].y*current_path[i].y;

        //remove last position from path
        current_path.pop_back();
        return;
    }

    for(int nextDir = 0; nextDir < 4; nextDir++){ //recursion step for all directions
        //dont go backwards
        if((dir + 2) % 4 != nextDir) saw_step(i + 1, nextDir);
    }
    //when done, remove step
    current_path.pop_back();
}
tuple saw(int steps)
{
    //reset global variables
    progress_counter = 0;
    Z = 0;
    endpoint_sum = 0;
    N_steps = steps;
    current_path.clear();
    current_path.emplace_back(0, 0);

    //1st recursion step
    for(int nextDir = 0; nextDir < 4; nextDir++){
        saw_step(1, nextDir);
    }

    current_path.pop_back(); //clear current path
    std::cout << "\e[2K \r" << std::flush; //clear progress line

    //calculate <R_ee^2>
    double mean_squared = (double) endpoint_sum / Z;


    return {Z, mean_squared};
}

int main(int argc, char *argv[])
{

    out << "# N\tZ\t<R_ee^2>" << std::endl;
    for(int N = 1; N < 21; N++){
        std::cout << "\nN = " << N << std::endl;

        tuple data = saw(N); //run self avoiding walk for given N

        std::cout << "Z = " << data.Z << std::endl;
        std::cout << "<R_ee^2> = " << data.mean_squared << std::endl;

        out << N << " " << data.Z << " " << data.mean_squared << std::endl;
    }

    out.close();
    return 0;
}