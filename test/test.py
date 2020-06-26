#!/usr/bin/env python3

import chemutils as cu


# TODO TODO probably try to do most tests with and without cache

VERBOSE = False

def test_ignore_cache():
    # TODO currently only aiming to test the 'if_null' case, but expand to cover
    # all values of ignore_cache
    # TODO manually edit temporary cache (maybe just disable all atexit writing
    # at start of these tests?) to guarantee it is in cache pointing to null
    null_cached_name = 'methyl (Z)-3-hexenoate'
    import ipdb; ipdb.set_trace()
    assert null_cached_name in cu.cache['name']['inchi'], \
        'test pre-condition failed'

    # TODO TODO test both case where non-normed points to null as well as 
    # normed points to null (or maybe do that in below, since it should only
    # matter in the try_non_normalized case?)
    import ipdb; ipdb.set_trace()


def test_try_non_normalized():
    # TODO could manually set the cache value to null or something maybe (for
    # the original, non-normalized key, in a case where normalization has an
    # effect, to test this independent of future improvements to the
    # normalization step (may need to do something else...)

    # Matt pointed out that this lookup fails, even though the PubChem REST
    # calls with the original should work. It seems to be because
    # `normalize_name` converts it to 'methyl(z)-3-hexenoate'
    orig_name = 'methyl (Z)-3-hexenoate'
    inchi = cu.convert(orig_name, from_type='name',
        try_non_normalized=True, verbose=VERBOSE, #ignore_cache=True
    )

    # TODO also test a case w/ it False

    import ipdb; ipdb.set_trace()


if __name__ == '__main__':
    VERBOSE = True
    # TODO TODO TODO really want this? just ignore_cache=True? clear_cache?
    #cu.delete_cache()
    # TODO maybe cu.clear_cache() here? (checking that it actually won't
    # autosave overwriting existing atexit...)
    test_try_non_normalized()
    #test_ignore_cache()

