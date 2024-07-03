#ifndef __EnumHash_h__
#define __EnumHash_h__

// Hash for Enum Class
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

template <typename Key, typename T>
using COMPASUnorderedMap = std::unordered_map<Key, T, HashType<Key>>;


#endif // __EnumHash_h__
