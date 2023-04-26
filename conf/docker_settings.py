"""
Django settings for iSkyLIMS project.

Generated by 'django-admin startproject' using Django 1.11.4.

For more information on this file, see
https://docs.djangoproject.com/en/1.11/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.11/ref/settings/
"""

import os

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = PLACEHOLDER

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

ALLOWED_HOSTS = ['localhost','127.0.0.1']

# Application definition

INSTALLED_APPS = [
    'core',
    'clinic',
    'iSkyLIMS_wetlab',
    'drylab',
    'django_utils',
    'mptt',
    'crispy_forms',
    'django_crontab',
    'django_mptt_admin',
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_extensions',
    'rest_framework',
    'django_cleanup', # should go after your apps
]

APPS_NAMES = [ ['clinic', 'Clinic'],
    ['iSkyLIMS_wetlab', 'Masive Sequencing'],
    ['drylab','Requesting Services']
    ]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'iSkyLIMS.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [BASE_DIR + '/documents/drylab/services_templates'],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
                'django.template.context_processors.i18n',
            ],
        },
    },
]

WSGI_APPLICATION = 'iSkyLIMS.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.11/ref/settings/#databases
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'USER': 'django',
        'PASSWORD':'djangopass',
        'PORT':'3306',
        'NAME': 'iSkyLIMS',
        'HOST':'db1',
        'TEST': {
            'NAME': 'iSkyLIMS_test',
        },
    },
}


# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/1.11/topics/i18n/

LANGUAGE_CODE = 'en-us'
#TIME_ZONE = 'UTC' --> 'Europe/Madrid' for forms such as 'blog entries'
TIME_ZONE = 'Europe/Madrid'

USE_I18N = True

USE_L10N = True

USE_TZ = False


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.11/howto/static-files/

STATIC_URL = '/static/'
STATIC_ROOT = os.path.join(BASE_DIR, "static/")

#  Media settings
MEDIA_URL = '/documents/'
MEDIA_ROOT = os.path.join(BASE_DIR, 'documents/')

## Crispy forms settings
CRISPY_TEMPLATE_PACK = "bootstrap3"

# Redirect to home URL after login (Default redirects to /accounts/profile/)
LOGIN_REDIRECT_URL = '/'

EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'  # During development only
#EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'

# EMAIL settings
EMAIL_HOST = 'localhost'
EMAIL_PORT = '25'
EMAIL_HOST_USER = 'bioinfo'
EMAIL_HOST_PASSWORD = ''
EMAIL_USE_TLS = False
EMAIL_ISKYLIMS = "iSkyLIMS@localhost"
ALLOWED_EMAIL_DOMAINS = []

LOG_CRONTAB_FILE = os.path.join(BASE_DIR, 'logs', 'crontab.log')
LOG_CLEAN_FILE = os.path.join(BASE_DIR, 'logs', 'crontab_cleanup.log')


# Crontab settings
CRONJOBS = [
    #('2-59/5 * * * *', 'iSkyLIMS_wetlab.cron.update_run_in_recorded_state',  '>>' + LOG_CRONTAB_FILE), # run every 5 min wit an offset of 2 minutes
    #('0 0 2 * *', 'iSkyLIMS_wetlab.cron.update_run_in_recorded_state',  '>>' + LOG_CRONTAB_FILE), # At minute 0 past every 2nd hour
    #('*/5 * * * *', 'iSkyLIMS_wetlab.cron.check_not_finish_run', '>>' + LOG_CRONTAB_FILE) # run every 5 min

    #('0 */2 * * *', 'iSkyLIMS_wetlab.cron.look_for_miseq_runs', '>>' + LOG_CRONTAB_FILE) # At minute 0 past every 2nd hour
    ('0 0 2 * *', 'iSkyLIMS_wetlab.cron.looking_for_new_runs', '>>' + LOG_CRONTAB_FILE),
    ('0 0 1 * *', 'iSkyLIMS_wetlab.cron.delete_invalid_run', '>>' + LOG_CLEAN_FILE)
    ]

CRONTAB_COMMAND_SUFFIX = '2>&1'

#SITE_ID =1
DEFAULT_AUTO_FIELD = 'django.db.models.AutoField'

