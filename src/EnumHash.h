#ifndef __EnumHash_h__
#define __EnumHash_h__

// functor object to calculate hash for enum class
// note: this was written because c++11 didn't allow the use of enum hashes,
// but I believe it was considered a defect and has been fixed in c++14
// check http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-defects.html#2148
// we may be able to remove this
struct EnumClassHash
{
    template <typename T>
    std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

template <typename Key>
using HashType = typename std::conditional<std::is_enum<Key>::value, EnumClassHash, std::hash<Key>>::type;

// type alias for unordered map with enum hash
template <typename Key, typename T>
using COMPASUnorderedMap = std::unordered_map<Key, T, HashType<Key>>;


#endif // __EnumHash_h__
