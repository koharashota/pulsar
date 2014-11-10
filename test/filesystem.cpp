#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace fs = boost::filesystem;

int main()
{
    try {
        const fs::path path("dir1/a.txt");
        const std::time_t last_update = fs::last_write_time(path);

        const boost::posix_time::ptime time = boost::posix_time::from_time_t(last_update);
        std::cout << time << std::endl;
    }
    catch (fs::filesystem_error& ex) {
        std::cout << "エラー発生！ : " << ex.what() << std::endl;
    }
}
